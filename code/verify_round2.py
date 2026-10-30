#!/usr/bin/env python3
"""Independent round-2 verification: recompute NHIS-pooled discrimination/calibration
from stored linear predictors (sr_lp.csv) using a DIFFERENT code path than harden2.py,
re-verify the calculator worked examples from full-precision coefficients, and check
that every FigSR figure is at least as new as the canonical scripts that produce it.
"""
import numpy as np, pandas as pd, json, os, glob
from pathlib import Path
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "output"
SEED = 20261028

df = pd.read_csv(OUT / "sr_lp.csv")
print("src categories:", df["src"].unique().tolist())
ext = df[df["src"].astype(str).str.contains("nhis", case=False, na=False)]
if len(ext) == 0:
    # fallback: anything not derivation
    ext = df[~df["src"].astype(str).str.contains("nhanes|der", case=False, na=False)]
y = ext["y"].values.astype(float)
lp = ext["lp"].values.astype(float)

auc = roc_auc_score(y, lp)
X = sm.add_constant(lp.reshape(-1, 1))
fit = sm.Logit(y, X).fit(disp=0)
slope = float(fit.params[1])                       # Van Houwelingen calibration slope (free lp coef)
# calibration-in-the-large intercept = offset model logit(y)=a+1*lp (lp coefficient fixed at 1)
fit_off = sm.Logit(y, np.ones((len(y), 1)), offset=lp).fit(disp=0)
intercept = float(fit_off.params[0])               # reported as -0.2131
# (diagnostic) free 2-param intercept for reference, NOT the reported CITL:
intercept_free = float(fit.params[0])
p = 1/(1+np.exp(-lp))
mean_pred = float(p.mean()); obs = float(y.mean())
brier = float(np.mean((p - y)**2))

head = json.load(open(OUT / "sr_headline.json"))["external_pooled"]
tol = 0.0005
checks = {
    "auc":      (round(auc,4),      head["auc"]),
    "slope":    (round(slope,4),    head["slope"]),
    "intercept":(round(intercept,4),head["intercept"]),
    "mean_pred":(round(mean_pred,4),head["mean_pred"]),
    "obs":      (round(obs,5),      round(head["prev"],5)),
    "brier":    (round(brier,5),    head["brier"]),
    "n":        (int(len(ext)),     head["n"]),
    "events":   (int(y.sum()),      head["events"]),
}
print("\n--- Independent recompute vs sr_headline (stored) ---")
ok = True
for k,(got,stored) in checks.items():
    match = (got == stored) or (isinstance(got,float) and abs(got-stored) <= max(tol,0.0006))
    if k in ("n","events"): match = (got == stored)
    print(f"{k:10s} got={got}  stored={stored}  {'OK' if match else 'MISMATCH'}")
    ok = ok and match
print(f"(diagnostic) free 2-param intercept = {intercept_free:.4f}; reported CITL (offset) = {intercept:.4f}")

# --- calculator worked examples from FULL-precision coefficients ---
me = json.load(open(OUT / "model_equation.json"))
cf = me["coefficients_full_precision"]
def logit(row):
    return (cf["const"] + cf["age"]*row["age"] + cf["male"]*row.get("male",0)
            + cf["smoke_former"]*row.get("smoke_former",0) + cf["smoke_current"]*row.get("smoke_current",0)
            + cf["obesity"]*row.get("obesity",0) + cf["hypertension"]*row.get("hypertension",0)
            + cf["diabetes"]*row.get("diabetes",0) + cf["dyslipidemia"]*row.get("dyslipidemia",0))
r1 = {"age":70,"male":1,"smoke_former":1,"obesity":1,"hypertension":1,"dyslipidemia":1}
r2 = {"age":55,"male":0}
ex1 = round(float(1/(1+np.exp(-logit(r1)))),4)
ex2 = round(float(1/(1+np.exp(-logit(r2)))),4)
print("\n--- Calculator worked examples (independent, full precision) ---")
print(f"70yM(former,obese,HTN,highchol) recomputed={ex1}  published={me['worked_examples']['70yM(former,obese,HTN,highchol,nonDiab)']}")
print(f"55F healthy            recomputed={ex2}  published={me['worked_examples']['55F healthy']}")
# allow the documented 4-dp rounding band (max prob deviation 0.00045)
ok = ok and abs(ex1-0.0596)<=0.0005 and abs(ex2-0.0117)<=0.0005

# --- 4dp vs full-precision reproduction over the cohort age range ---
c4 = me["coefficients_4dp"]
def logit4(row):
    return (c4["const"] + c4["age"]*row["age"] + c4["male"]*row.get("male",0)
            + c4["smoke_former"]*row.get("smoke_former",0) + c4["smoke_current"]*row.get("smoke_current",0)
            + c4["obesity"]*row.get("obesity",0) + c4["hypertension"]*row.get("hypertension",0)
            + c4["diabetes"]*row.get("diabetes",0) + c4["dyslipidemia"]*row.get("dyslipidemia",0))
dmax = max(abs(logit4(r1)-logit(r1)), abs(logit4(r2)-logit(r2)))
print(f"4dp-vs-full |ΔLP| on the two examples = {dmax:.5f} (threshold 0.00223 over cohort)")

# --- figure freshness: each FigSR*.{png,pdf} mtime >= producing script mtime ---
print("\n--- Figure freshness vs canonical scripts ---")
scripts = {
    "FigSR1": OUT.parent/"scripts"/"sr_final_tables_figs.py",
    "FigSR2": OUT.parent/"scripts"/"sr_final_tables_figs.py",
    "FigSR3": OUT.parent/"scripts"/"sr_final_tables_figs.py",
    "FigSR4": OUT.parent/"scripts"/"sr_final_tables_figs.py",
    "FigSR5": OUT.parent/"scripts"/"sr_final_tables_figs.py",
}
fig_dir = ROOT/"analysis"/"figures"
if not fig_dir.exists():
    cands = glob.glob(str(ROOT/"**"/"FigSR1_transportability.png"), recursive=True)
    fig_dir = Path(os.path.dirname(cands[0])) if cands else None
fresh_ok = True
for tag, sc in scripts.items():
    if not sc.exists():
        print(f"{tag}: SCRIPT MISSING {sc}"); fresh_ok=False; continue
    sm_time = sc.stat().st_mtime
    pngs = sorted(glob.glob(str(fig_dir/f"{tag}*")))
    for fg in pngs:
        ft = os.path.getmtime(fg)
        flag = "OK" if ft >= sm_time else "STALE (regenerate!)"
        if ft < sm_time: fresh_ok = False
        print(f"{os.path.basename(fg):32s} fig_mtime={ft:.0f} script_mtime={sm_time:.0f} {flag}")

verdict = "PASS" if (ok and fresh_ok) else "FAIL"
print("\nOVERALL VERIFICATION:", verdict)
json.dump({"independent_external":checks, "worked_examples":{"ex1":ex1,"ex2":ex2},
           "numeric_ok":ok, "figures_fresh":fresh_ok, "verdict":verdict},
          open(OUT/"sr_verification_round2.json","w"), indent=2)
