"""ROUND-3 independent verification of the design-weighted external calibration.
Re-derives pooled weighted AUROC / slope / CITL / mean-pred / obs from nhis_external.csv
using an INDEPENDENT weighted-AUC implementation (Skoog rank trick) and a plain weighted
logistic via IRLS-free statsmodels GLM, then compares to sr_weighted_calibration.json.
Also checks FigSR6 freshness vs its producing script.
"""
import numpy as np, pandas as pd, json, os, glob
from pathlib import Path
import statsmodels.api as sm
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"; FIG=ROOT/"figures"

ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
ORDER=["age","male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
COEF=json.load(open(OUT/"model_equation.json"))["coefficients_full_precision"]
B=np.array([COEF["const"]]+[COEF[k] for k in ORDER])
lp=np.column_stack([np.ones(len(ext)),ext[ORDER].values.astype(float)])@B
y=ext["amd_self"].values.astype(float); w=ext["weight"].values.astype(float); p=1/(1+np.exp(-lp))

def weighted_auc(y,s,w):
    # Skoog (2008) weighted AUC via weighted Mann-Whitney
    pos=y==1; neg=~pos
    wp=w[pos]; wn=w[neg]; sp=s[pos]; sn=s[neg]
    # for each positive, weight of negatives strictly below + 0.5*ties
    sn_sorted=np.sort(sn); wn_cum=np.cumsum(wn[np.argsort(sn)])/1.0
    order=np.argsort(sn); sn_s=sn[order]; wn_s=wn[order]
    below=[]; tie=[]
    for x in sp:
        lo=sn_s< x; eq=sn_s==x
        below.append(wn_s[lo].sum()); tie.append(wn_s[eq].sum())
    num=np.sum(wp*(np.array(below)+0.5*np.array(tie)))
    den=np.sum(wp)*np.sum(wn)
    return num/den

wauc=weighted_auc(y,lp,w)
wslope=float(sm.GLM(y,sm.add_constant(lp),family=sm.families.Binomial(),freq_weights=w).fit().params[1])
wcitl=float(sm.GLM(y,np.ones((len(y),1)),family=sm.families.Binomial(),offset=lp,freq_weights=w).fit().params[0])
wmp=float(np.average(p,weights=w)); wob=float(np.average(y,weights=w))

stored=json.load(open(OUT/"sr_weighted_calibration.json"))["pooled"]["weighted"]
checks={"w_auc":(round(wauc,4),stored["auc"]),"w_slope":(round(wslope,4),stored["slope"]),
        "w_citl":(round(wcitl,4),stored["citl"]),"w_meanpred":(round(wmp,4),stored["mean_pred"]),
        "w_obs":(round(wob,4),stored["obs"]),"w_overpred":(round(wmp/wob,3),stored["overpred_ratio"])}
ok=True
print("--- Independent weighted recompute vs stored ---")
for k,(g,s) in checks.items():
    match=abs(g-s)<=0.0006; ok=ok and match
    print(f"{k:12s} got={g}  stored={s}  {'OK' if match else 'MISMATCH'}")

# figure freshness for FigSR6 vs producing script
sc=ROOT/"analysis"/"scripts"/"sr_weighted_calibration.py"
sm_time=sc.stat().st_mtime; fresh=True
for fg in sorted(glob.glob(str(FIG/"FigSR6*"))):
    ft=os.path.getmtime(fg); flag="OK" if ft>=sm_time else "STALE"; fresh=fresh and (ft>=sm_time)
    print(f"Fig6 {os.path.basename(fg):34s} {'OK' if ft>=sm_time else 'STALE'}")
verdict="PASS" if (ok and fresh) else "FAIL"
print("\nROUND-3 WEIGHTED VERIFICATION:",verdict)
json.dump({"independent_weighted":checks,"numeric_ok":ok,"fig6_fresh":fresh,"verdict":verdict},
          open(OUT/"sr_verification_round3.json","w"),indent=2)
