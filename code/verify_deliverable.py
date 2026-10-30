"""Pre-delivery verification (independent recompute + traceability + honesty grep)."""
import os, json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
OUT=r"H:\眼科论文\题目4\analysis\output"; MS=r"H:\眼科论文\题目4\manuscript"
import re
find=json.load(open(os.path.join(OUT,"findings.json")))
extlp=pd.read_csv(os.path.join(OUT,"external_lp.csv"))
pp=extlp[extlp.cohort=="NHIS_pooled"]
y=pp['y'].values; lp=pp['lp'].values
auc_ind=roc_auc_score(y,lp)
slope_ind=sm.Logit(y.astype(float),sm.add_constant(lp)).fit(disp=0).params[1]
f=find['external'][0]
print("Independent NHIS-pooled AUROC:",round(auc_ind,4),"| findings:",f['auc'])
print("Independent calibration slope:",round(slope_ind,4),"| findings:",f['calib_slope'])
assert abs(auc_ind-f['auc'])<1e-3, "AUROC mismatch"
assert abs(slope_ind-f['calib_slope'])<1e-2, "slope mismatch"
# cross-check Tables.md contains the same headline numbers as findings
tables=open(os.path.join(OUT,"Tables.md"),encoding='utf-8').read()
for token in ["0.775","0.848","-0.233","0.782","0.749","0.894"]:
    assert token in tables or token in open(os.path.join(MS,"manuscript.md"),encoding='utf-8').read(), f"missing {token}"
print("Tables/manuscript contain headline numbers: OK")
# honesty grep: manuscript must not claim ELSA/SHARE/SAGE were validated
ms=open(os.path.join(MS,"manuscript.md"),encoding='utf-8').read()
for w in ["ELSA","SHARE","SAGE","BRFSS"]:
    idx=ms.find(w)
    print(f"  '{w}' present:", idx>=0, "(only in scope/limitations context expected)")
assert "validated in ELSA" not in ms and "validated in SHARE" not in ms
# traceability: figures newer than analysis scripts that produce their data
figs=sorted(os.listdir(r"H:\眼科论文\题目4\figures"))
mtime_max_script=max(os.path.getmtime(os.path.join(r"H:\眼科论文\题目4\analysis\scripts",f)) for f in os.listdir(r"H:\眼科论文\题目4\analysis\scripts") if f.endswith(".py"))
older=[f for f in figs if os.path.getmtime(os.path.join(r"H:\眼科论文\题目4\figures",f))<mtime_max_script]
print("Figures older than newest script (STALE):",older if older else "none")
print("n figure files:",len([f for f in figs if f.endswith(('.png','.pdf'))]))
print("\nVERIFICATION PASSED")
