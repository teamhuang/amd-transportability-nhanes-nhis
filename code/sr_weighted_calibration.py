"""ROUND-3 core addition: DESIGN-WEIGHTED external calibration.
Recompute, under the NHIS complex-survey weights, every headline calibration quantity
(AUROC, Van Houwelingen slope, calibration-in-the-large, mean predicted vs observed,
Brier) plus a weighted cross-year intercept update and weighted decision curve, for the
frozen full model. Reported SIDE BY SIDE with the unweighted (primary) estimand so the
manuscript can bound the drift/repair claim to the target population and estimator.
Frozen coefficients are read from model_equation.json (full precision) so nothing is refit.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"; FIG=ROOT/"figures"
import numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, brier_score_loss
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ORDER=["age","male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
COEF=json.load(open(OUT/"model_equation.json"))["coefficients_full_precision"]
BVEC=np.array([COEF["const"]]+[COEF[k] for k in ORDER])   # const + ordered predictors

ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
def lp_of(df):
    M=np.column_stack([np.ones(len(df)),df[ORDER].values.astype(float)])
    return M@BVEC
def probs(lp): return 1/(1+np.exp(-lp))

def calib(y,lp,w=None):
    """w=None -> unweighted (matches canonical); else design-weighted point estimates."""
    p=probs(lp); y=np.asarray(y,float); lp=np.asarray(lp,float)
    sw=None if w is None else np.asarray(w,float)
    auc=float(roc_auc_score(y,lp,sample_weight=sw))
    X=sm.add_constant(lp)
    if w is None:
        slope=float(sm.Logit(y,X).fit(disp=0).params[1])
        citl=float(sm.Logit(y,np.ones((len(y),1)),offset=lp).fit(disp=0).params[0])
        mean_pred=float(p.mean()); obs=float(y.mean()); brier=float(brier_score_loss(y,p))
    else:
        slope=float(sm.GLM(y,X,family=sm.families.Binomial(),freq_weights=sw).fit().params[1])
        citl=float(sm.GLM(y,np.ones((len(y),1)),family=sm.families.Binomial(),
                          offset=lp,freq_weights=sw).fit().params[0])
        mean_pred=float(np.average(p,weights=sw)); obs=float(np.average(y,weights=sw))
        brier=float(np.average((p-y)**2,weights=sw))
    return dict(auc=round(auc,4),slope=round(slope,4),citl=round(citl,4),
                mean_pred=round(mean_pred,4),obs=round(obs,4),brier=round(brier,5),
                overpred_ratio=round(mean_pred/obs,3) if obs>0 else None)

lp=lp_of(ext); y=ext["amd_self"].values; w=ext["weight"].values
res={"pooled":{"unweighted":calib(y,lp),"weighted":calib(y,lp,w)}}
for yr in (2016,2017):
    m=ext["survey_year"].values==yr
    res[f"NHIS_{yr}"]={"unweighted":calib(y[m],lp[m]),"weighted":calib(y[m],lp[m],w[m])}

# ---- cross-year intercept update, under each estimand ----
def intercept_on(idx):
    """calibration-in-the-large intercept on the selected rows (idx bool mask)."""
    if idx is None:
        return float(sm.Logit(y,np.ones((len(y),1)),offset=lp).fit(disp=0).params[0])
    return float(sm.Logit(y[idx],np.ones((int(idx.sum()),1)),offset=lp[idx]).fit(disp=0).params[0])
def intercept_on_w(idx):
    return float(sm.GLM(y[idx],np.ones((int(idx.sum()),1)),family=sm.families.Binomial(),
                        offset=lp[idx],freq_weights=w[idx]).fit().params[0])
def cross_update(estim):
    out={}
    for src,tgt in [(2016,2017),(2017,2016)]:
        ms=ext["survey_year"].values==src; mt=ext["survey_year"].values==tgt
        if estim=="unweighted":
            a=intercept_on(ms)
            upd=probs(lp[mt]+a)
            mp=float(upd.mean()); ob=float(y[mt].mean()); br=float(np.mean((upd-y[mt])**2))
        else:
            a=intercept_on_w(ms)
            upd=probs(lp[mt]+a)
            mp=float(np.average(upd,weights=w[mt])); ob=float(np.average(y[mt],weights=w[mt]))
            br=float(np.average((upd-y[mt])**2,weights=w[mt]))
        frozen=res[f"NHIS_{tgt}"][estim]
        out[f"update_{src}_to_{tgt}"]=dict(intercept=round(a,4),
            frozen_mean_pred=frozen["mean_pred"],updated_mean_pred=round(mp,4),
            observed=round(ob,4),frozen_brier=frozen["brier"],updated_brier=round(br,5))
    return out
res["cross_year_update"]={"unweighted":cross_update("unweighted"),
                          "weighted":cross_update("weighted")}

# ---- weighted vs unweighted decision curve (referral) at prespecified thresholds ----
def dca(y,p,w,th_list):
    y=np.asarray(y,float); p=np.asarray(p,float); out=[]
    W=np.ones_like(y) if w is None else np.asarray(w,float)
    nn=float(np.sum(W))
    prev=float(np.sum(W*y)/nn)                       # (weighted) prevalence
    for pt in th_list:
        flag=p>=pt
        tp=float(np.sum(W[(flag)&(y==1)])); fp=float(np.sum(W[(flag)&(y==0)]))
        nb=(tp-(fp*pt/(1-pt)))/nn
        nb_all=prev-(1-prev)*(pt/(1-pt))             # treat-all net benefit, correct form
        out.append(dict(threshold=pt,
                        flagged_per_1000=round(1000*float(np.sum(W[flag]))/nn,1),
                        cases_per_1000=round(1000*tp/nn,1),
                        nb=round(nb,5), nb_treat_all=round(nb_all,5)))
    return out
TH=[0.02,0.05,0.08,0.10]; pf=probs(lp)
res["dca"]={"unweighted":dca(y,pf,None,TH),"weighted":dca(y,pf,w,TH)}

json.dump(res,open(OUT/"sr_weighted_calibration.json","w"),indent=2)
print(json.dumps(res,indent=1))

# ---- comparison table + figure ----
rows=[]
for k in ["pooled","NHIS_2016","NHIS_2017"]:
    u=res[k]["unweighted"]; v=res[k]["weighted"]
    for q in ["auc","slope","citl","mean_pred","obs","overpred_ratio"]:
        rows.append(dict(set=k,quantity=q,unweighted=u[q],weighted=v[q],diff=round((v[q]-u[q]),4) if isinstance(u[q],(int,float)) and isinstance(v[q],(int,float)) else None))
tab=pd.DataFrame(rows); tab.to_csv(OUT/"sr_weighted_vs_unweighted.csv",index=False)

if os.environ.get("LEGACY_FIGURES", "0") != "1":
    print("[legacy figure block skipped; figures are produced by make_figures_nature_round17.py]")
    raise SystemExit(0)

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8})
OI=dict(blue="#0072B2",orange="#E69F00",green="#009E73",verm="#D55E00",purple="#CC79A7",grey="#7f7f7f",black="#000000")
f,(a1,a2)=plt.subplots(1,2,figsize=(6.6,3.2))
sets=["pooled","NHIS_2016","NHIS_2017"]; x=np.arange(3)
uw=[res[s]["unweighted"]["mean_pred"] for s in sets]; ww=[res[s]["weighted"]["mean_pred"] for s in sets]
uo=[res[s]["unweighted"]["obs"] for s in sets];      wo=[res[s]["weighted"]["obs"] for s in sets]
a1.bar(x-0.2,uw,0.38,color=OI["verm"],label="pred (unweighted)"); a1.bar(x-0.2+0.18,ww,0.38,color=OI["purple"],label="pred (weighted)")
a1.bar(x+0.2,uo,0.38,color=OI["blue"],label="obs (unweighted)"); a1.bar(x+0.2+0.18,wo,0.38,color=OI["green"],label="obs (weighted)")
a1.set_xticks(x);a1.set_xticklabels([s.replace("NHIS_","") for s in sets]);a1.set_ylabel("Risk (probability)")
a1.legend(frameon=False,fontsize=5.5);a1.set_title("A  Predicted vs observed, by estimand",fontsize=8)
us=[res[s]["unweighted"]["citl"] for s in sets]; ws=[res[s]["weighted"]["citl"] for s in sets]
a2.bar(x-0.2,us,0.38,color=OI["verm"],label="CITL (unweighted)"); a2.bar(x+0.2,ws,0.38,color=OI["purple"],label="CITL (weighted)")
a2.axhline(0,color=OI["black"],lw=0.6);a2.set_xticks(x);a2.set_xticklabels([s.replace("NHIS_","") for s in sets])
a2.set_ylabel("Calibration-in-the-large");a2.legend(frameon=False,fontsize=6);a2.set_title("B  Intercept depends on estimand",fontsize=8)
f.tight_layout()
f.savefig(FIG/"FigSR6_weighted_calibration.png",dpi=300,bbox_inches="tight")
f.savefig(FIG/"FigSR6_weighted_calibration.pdf",bbox_inches="tight"); plt.close(f)
print("SAVED sr_weighted_calibration.json, sr_weighted_vs_unweighted.csv, FigSR6")
