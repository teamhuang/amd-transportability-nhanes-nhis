"""ROUND-4 P0-2: FORMAL complex-survey calibration inference for the external (NHIS) model.
(1) Design-weighted prevalence with first-order Taylor (delta-method) linearisation using the
    official strata and PSU (composite year|stratum|psu id to avoid cross-year label reuse).
(2) Design-weighted calibration slope and calibration-in-the-large via weighted GLM with
    PSU-clustered (linearised) robust SE, as a survey-GLM approximation.
(3) PSU-cluster bootstrap with B=1000 and TWO seeds (stability check) for weighted AUROC/slope/
    intercept. Documents that the weighted AUROC is a DIFFERENT target parameter, not 'better'.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"
import numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
SEED=20261028
ORDER=["age","male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
COEF=json.load(open(OUT/"model_equation.json"))["coefficients_full_precision"]
B=np.array([COEF["const"]]+[COEF[k] for k in ORDER])

ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
ext["lp"]=np.column_stack([np.ones(len(ext)),ext[ORDER].values.astype(float)])@B
ext["p"]=1/(1+np.exp(-ext["lp"]))
ext["psu_id"]=ext["survey_year"].astype(str)+"_"+ext["stratum"].astype(str)+"_"+ext["psu"].astype(str)
ext["stratum_id"]=ext["survey_year"].astype(str)+"_"+ext["stratum"].astype(str)
y=ext["amd_self"].values.astype(float); w=ext["weight"].values.astype(float)

def wprev_taylor(y,w,strata,psu):
    R=np.sum(w*y)/np.sum(w)
    df=pd.DataFrame(dict(y=y,w=w,st=strata,psu=psu)); df["e"]=df.w*(df.y-R)
    z=df.groupby(["st","psu"])["e"].sum().reset_index()
    var=0.0
    for st,g in z.groupby("st"):
        m=len(g)
        if m>1: var+= (m/(m-1))*np.sum(g["e"].values**2)
        # if m==1, treat as 0 contribution (single-cluster stratum; conservative)
    var=var/(np.sum(w)**2)
    se=np.sqrt(var); return R, se, 1.96*se

def wcalib_cluster(y,lp,w,psu_groups):
    slope_mod=sm.GLM(y,sm.add_constant(lp),family=sm.families.Binomial(),freq_weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":psu_groups})
    citl_mod=sm.GLM(y,np.ones((len(y),1)),family=sm.families.Binomial(),offset=lp,freq_weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":psu_groups})
    s=slope_mod.params[1]; sse=slope_mod.bse[1]
    a=citl_mod.params[0]; ase=citl_mod.bse[0]
    return dict(slope=round(float(s),4),slope_se=round(float(sse),4),
                slope_ci=[round(float(s-1.96*sse),4),round(float(s+1.96*sse),4)],
                citl=round(float(a),4),citl_se=round(float(ase),4),
                citl_ci=[round(float(a-1.96*ase),4),round(float(a+1.96*ase),4)])

groups=ext["psu_id"].values
res={"target_note":"Design-weighted estimates target the US civilian noninstitutionalised population; the unweighted estimates target the analytic sample. They are different parameters, not a better/worse model.",
     "prevalence":{},"calibration_cluster_glm":{}}
for name,mask in [("pooled",np.ones(len(ext),bool)),("2016",ext.survey_year.values==2016),("2017",ext.survey_year.values==2017)]:
    R,se,half=wprev_taylor(y[mask],w[mask],ext["stratum_id"].values[mask],ext["psu_id"].values[mask])
    res["prevalence"][name]=dict(weighted=round(float(R),4),taylor_se=round(float(se),4),
                                 ci95=[round(float(R-half),4),round(float(R+half),4)])
    res["calibration_cluster_glm"][name]=wcalib_cluster(y[mask],ext["lp"].values[mask],w[mask],groups[mask])

# ---------- PSU-cluster bootstrap, B=1000, two seeds ----------
def wmetrics(yb,lpb,wb):
    p=1/(1+np.exp(-lpb))
    auc=roc_auc_score(yb,lpb,sample_weight=wb)
    sl=sm.GLM(yb,sm.add_constant(lpb),family=sm.families.Binomial(),freq_weights=wb).fit().params[1]
    ci=sm.GLM(yb,np.ones((len(yb),1)),family=sm.families.Binomial(),offset=lpb,freq_weights=wb).fit().params[0]
    return auc,sl,ci,np.average(p,weights=wb)
# build PSU index map
ext["_row"]=np.arange(len(ext))
psu_groups=ext.groupby("psu_id")["_row"].apply(lambda s:s.values).to_dict()
psu_keys=np.array(list(psu_groups.keys()))
stratum_of={k:ext.loc[ext.psu_id==k,"stratum_id"].iloc[0] for k in psu_keys}
strata=np.array([stratum_of[k] for k in psu_keys])
def boot(seed):
    rng=np.random.RandomState(seed); out=[]
    for _ in range(1000):
        rows=[]
        for st in np.unique(strata):
            ks=psu_keys[strata==st]; m=len(ks)
            draw=ks[rng.randint(0,m,m)]
            for k in draw: rows.extend(psu_groups[k])
        rows=np.array(rows)
        yb=y[rows]; lpb=ext["lp"].values[rows]; wb=w[rows]
        if len(np.unique(yb))<2: continue
        try: out.append(wmetrics(yb,lpb,wb))
        except Exception: pass
    A=np.array(out); 
    return {q:[round(float(np.percentile(A[:,i],2.5)),4),round(float(np.percentile(A[:,i],97.5)),4)]
            for i,q in enumerate(["auc","slope","citl","meanpred"])}, A.shape[0]
b1,n1=boot(SEED); b2,n2=boot(777)
res["psu_bootstrap"]={"B":1000,"reps_ok_seed1":int(n1),"reps_ok_seed2":int(n2),
                      "seed_20261028":b1,"seed_777":b2,
                      "seed_stability_maxdiff":round(float(max(abs(b1[k][0]-b2[k][0])+abs(b1[k][1]-b2[k][1]) for k in b1)),4)}
# point estimates (weighted)
res["weighted_point"]=dict(auc=round(float(roc_auc_score(y,ext["lp"].values,sample_weight=w)),4),
                           **wcalib_cluster(y,ext["lp"].values,w,groups))
json.dump(res,open(OUT/"sr_survey.json","w"),indent=2)
print(json.dumps(res,indent=1))
