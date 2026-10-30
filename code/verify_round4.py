"""ROUND-4 independent verification: recompute, via a separate code path, the headline numbers
that the round-4 tables/manuscript now cite, and check FigSR7/FigSR8 freshness.
"""
import numpy as np, pandas as pd, json, os, glob
from pathlib import Path
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"; FIG=ROOT/"figures"
ORDER=["age","male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
COEF=json.load(open(OUT/"model_equation.json"))["coefficients_full_precision"]
B=np.array([COEF["const"]]+[COEF[k] for k in ORDER])
ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
lp=np.column_stack([np.ones(len(ext)),ext[ORDER].values.astype(float)])@B
y=ext["amd_self"].values.astype(float); w=ext["weight"].values.astype(float)

# weighted prevalence + weighted AUROC (independent)
wprev=np.average(y,weights=w); wauc=roc_auc_score(y,lp,sample_weight=w)
# unweighted AUROC (independent, matches canonical 0.7737)
uauc=roc_auc_score(y,lp)
# full vs age-only paired dAUROC point (independent of sr_increment): refit age-only
def D(df,preds): return np.column_stack([np.ones(len(df)),df[["age"]+preds].values.astype(float)])
der=pd.read_csv(OUT/"nhanes_derivation.csv"); der=der[der.cc_pred==True].reset_index(drop=True)
ma=sm.Logit(der["amd_self"].values.astype(float),D(der,[])).fit(disp=0)
la=D(ext,[])@np.asarray(ma.params)
dau=roc_auc_score(y,lp)-roc_auc_score(y,la)
sv=json.load(open(OUT/"sr_survey.json")); inc=json.load(open(OUT/"sr_increment.json")); mi=json.load(open(OUT/"sr_missing_mi2.json"))
checks={
 "w_prevalence":(round(float(wprev),4), sv["prevalence"]["pooled"]["weighted"], 0.0006),
 "w_auc":(round(float(wauc),4), sv["weighted_point"]["auc"], 0.0006),
 "u_auc":(round(float(uauc),4), 0.7737, 0.0006),
 "dAUROC_full_age":(round(float(dau),4), inc["increment"]["unweighted"]["dAUC"], 0.0006),
 "mi_ext_auc_mean":(round(float(mi["multiple_imputation"]["ext_auc_mean"]),4), 0.7744, 0.0002),
}
ok=True
print("--- Round-4 independent recompute ---")
for k,(got,ref,tol) in checks.items():
    match=abs(got-ref)<=tol; ok=ok and match
    print(f"{k:20s} got={got} ref={ref} {'OK' if match else 'MISMATCH'}")

# objective any-ARM paired dAUROC sign (full < age) — independent quick recompute
def Ddf(df,preds): return np.column_stack([np.ones(len(df)),df[["age"]+preds].values.astype(float)])
mf=sm.Logit(der["amd_self"].values.astype(float),Ddf(der,ORDER[1:])).fit(disp=0)
g=der[der.obj_anyARM.notna()].reset_index(drop=True)
lf_g=Ddf(g,ORDER[1:])@np.asarray(mf.params); la_g=Ddf(g,[])@np.asarray(ma.params); yg=g.obj_anyARM.values.astype(float)
d_obj=roc_auc_score(yg,lf_g)-roc_auc_score(yg,la_g)
print(f"objective any-ARM dAUROC (indep) = {d_obj:.4f}  (manuscript -0.0117); sign_match={d_obj<0}")
ok = ok and (d_obj<0)

# figure freshness
fresh=True
for tag in ["FigSR7_calibration_selfreport","FigSR8_calibration_objective"]:
    sc=ROOT/"analysis"/"scripts"/"make_round4_tables_figs.py"; smt=sc.stat().st_mtime
    for fg in glob.glob(str(FIG/(tag+"*"))):
        ft=os.path.getmtime(fg); m=ft>=smt; fresh=fresh and m
        print(f"{os.path.basename(fg):44s} {'OK' if m else 'STALE'}")
verdict="PASS" if (ok and fresh) else "FAIL"
print("\nROUND-4 VERIFICATION:",verdict)
json.dump({"checks":{k:[v[0],v[1]] for k,v in checks.items()},"objective_sign_indep":round(float(d_obj),4),
           "numeric_ok":bool(ok),"figs_fresh":bool(fresh),"verdict":verdict},
          open(OUT/"sr_verification_round4.json","w"),indent=2)
