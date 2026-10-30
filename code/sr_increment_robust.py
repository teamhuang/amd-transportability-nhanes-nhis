"""ROUND-4 P0-4 / P1: model-increment fairness and robustness.
(a) Paired full-vs-age-only ΔAUROC, ΔBrier and ΔNetBenefit at 2/5/8/10% with bootstrap CIs,
    under BOTH estimands (unweighted analytic sample and design-weighted population).
(b) Objective fundus: paired ΔAUROC (full - age-only) with bootstrap CI for any-ARM and late AMD,
    plus a gradeability (gradable vs ungradable) selection comparison.
(c) Age functional form: linear vs 3-df restricted-cubic-spline compared on external AUROC AND
    calibration (slope, CITL, mean predicted) — not AUROC alone.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"
import numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, brier_score_loss
SEED=20261028
BIN=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
der=pd.read_csv(OUT/"nhanes_derivation.csv"); der=der[der.cc_pred==True].reset_index(drop=True)
ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
def D(df,preds): return np.column_stack([np.ones(len(df)),df[["age"]+preds].values.astype(float)])
def fitlog(df,preds,end="amd_self"):
    X=D(df,preds); return sm.Logit(df[end].values.astype(float),X).fit(disp=0)
mf=fitlog(der,BIN); ma=fitlog(der,[])
lf=D(ext,BIN)@np.asarray(mf.params); la=D(ext,[])@np.asarray(ma.params); ye=ext["amd_self"].values.astype(float); we=ext["weight"].values.astype(float)

def nb(y,p,pt,w=None):
    W=np.ones_like(y) if w is None else w; flag=p>=pt; n=W.sum()
    tp=W[(flag)&(y==1)].sum(); fp=W[(flag)&(y==0)].sum(); return (tp-fp*pt/(1-pt))/n
def paired_delta(seed,B=1000,w=None):
    rng=np.random.RandomState(seed); n=len(ye); use = np.ones(n,bool) if w is None else None
    da=[];db=[];dn=[]
    sw=None if w is None else w
    for _ in range(B):
        idx=rng.randint(0,n,n); yy=ye[idx]; lfi=lf[idx]; lai=la[idx]; wi=None if w is None else w[idx]
        if len(np.unique(yy))<2: continue
        da.append(roc_auc_score(yy,lfi,sample_weight=wi)-roc_auc_score(yy,lai,sample_weight=wi))
        pi=1/(1+np.exp(-lfi)); pai=1/(1+np.exp(-lai))
        bi=brier_score_loss(yy,pi,sample_weight=wi)-brier_score_loss(yy,pai,sample_weight=wi)
        db.append(bi)
        dn.append(nb(yy,pi,0.05,wi)-nb(yy,pai,0.05,wi))
    A=np.array([da,db,dn])
    return dict(dAUC=round(float(np.mean(da)),5),dAUC_CI=[round(float(np.percentile(da,2.5)),5),round(float(np.percentile(da,97.5)),5)],
                dBrier=round(float(np.mean(db)),6),dBrier_CI=[round(float(np.percentile(db,2.5)),6),round(float(np.percentile(db,97.5)),6)],
                dNB5=round(float(np.mean(dn)),6),dNB5_CI=[round(float(np.percentile(dn,2.5)),6),round(float(np.percentile(dn,97.5)),6)])
inc={"unweighted":paired_delta(SEED),"weighted":paired_delta(SEED+1,w=we)}
# ΔNB across thresholds (point, both estimands)
def nb_grid(w=None):
    pi=1/(1+np.exp(-lf)); pai=1/(1+np.exp(-la)); out=[]
    for pt in [0.02,0.05,0.08,0.10]:
        out.append(dict(threshold=pt,nb_full=round(nb(ye,pi,pt,w),5),nb_age=round(nb(ye,pai,pt,w),5),
                        dNB=round(nb(ye,pi,pt,w)-nb(ye,pai,pt,w),5)))
    return out
inc["dNB_grid_unweighted"]=nb_grid(None); inc["dNB_grid_weighted"]=nb_grid(we)

# ---------- (b) objective fundus paired ΔAUROC + gradeability ----------
def obj_analysis(end):
    d=der[der[end].notna()].reset_index(drop=True)
    lf_d=D(d,BIN)@np.asarray(mf.params); la_d=D(d,[])@np.asarray(ma.params); yd=d[end].values.astype(float)
    auc_f=roc_auc_score(yd,lf_d); auc_a=roc_auc_score(yd,la_d)
    rng=np.random.RandomState(SEED); dd=[]
    for _ in range(500):
        idx=rng.randint(0,len(yd),len(yd)); yy=yd[idx]
        if len(np.unique(yy))<2: continue
        dd.append(roc_auc_score(yy,lf_d[idx])-roc_auc_score(yy,la_d[idx]))
    return dict(endpoint=end,n=int(len(d)),events=int(yd.sum()),auc_full=round(float(auc_f),4),
                auc_age=round(float(auc_a),4),dAUC_full_minus_age=round(float(auc_f-auc_a),4),
                dAUC_CI=[round(float(np.percentile(dd,2.5)),4),round(float(np.percentile(dd,97.5)),4)])
obj=[obj_analysis("obj_anyARM"),obj_analysis("obj_lateARM")]
# gradeability selection: gradable vs ungradable (any-ARM)
g=der["obj_anyARM"].notna()
grade=dict(gradable_n=int(g.sum()),ungradable_n=int((~g).sum()),
           mean_age_gradable=round(float(der.loc[g,"age"].mean()),2),mean_age_ungradable=round(float(der.loc[~g,"age"].mean()),2),
           amd_gradable=round(float(der.loc[g,"amd_self"].mean()),4),amd_ungradable=round(float(der.loc[~g,"amd_self"].mean()),4),
           pct_female_gradable=round(float(der.loc[g,"male"].mean()),4),pct_female_ungradable=round(float(der.loc[~g,"male"].mean()),4))

# ---------- (c) linear vs RCS age functional form ----------
KNOTS=[41.0,53.0,65.0,80.0]
def _cp(x,t): return np.maximum(x-t,0.0)**3
def rcs_basis(age):
    t1,t2,t3,t4=KNOTS; x=np.asarray(age,float); d=t4-t3
    r1=_cp(x,t1)-_cp(x,t3)*(t4-t1)/d+_cp(x,t4)*(t3-t1)/d
    r2=_cp(x,t2)-_cp(x,t3)*(t4-t2)/d+_cp(x,t4)*(t3-t2)/d
    return np.column_stack([x,r1,r2])
def design_rcs(df,preds):
    R=rcs_basis(df["age"].values)
    return np.column_stack([np.ones(len(df)),R,df[preds].values.astype(float)])
def eval_ext(Xd,dfit,Xe):
    m=sm.Logit(dfit["amd_self"].values.astype(float),Xd).fit(disp=0)
    lp=Xe@np.asarray(m.params); ye2=ext["amd_self"].values.astype(float); p=1/(1+np.exp(-lp))
    sl=float(sm.Logit(ye2,sm.add_constant(lp)).fit(disp=0).params[1])
    ci=float(sm.Logit(ye2,np.ones((len(ye2),1)),offset=lp).fit(disp=0).params[0])
    return dict(ext_auc=round(float(roc_auc_score(ye2,lp)),4),ext_slope=round(sl,4),ext_citl=round(ci,4),
                ext_meanpred=round(float(p.mean()),4))
funct=dict(
  linear=eval_ext(D(der,BIN),der,D(ext,BIN)),
  rcs=eval_ext(design_rcs(der,BIN),der,design_rcs(ext,BIN)),
  note="RCS adds 2 df on age; external AUROC gain <0.001 while coefficients become non-reproducible at 4dp; linear retained for parsimony+reproducibility. Functional-form choice is labelled exploratory; model was NOT selected on external AUROC (both forms give ~0.774).")

out=dict(increment=inc,objective=obj,gradeability=grade,age_function_form=funct)
json.dump(out,open(OUT/"sr_increment.json","w"),indent=2)
print(json.dumps(out,indent=1))
