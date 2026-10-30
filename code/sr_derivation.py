"""Derivation-layer table + internal validation for the LINEAR-age canonical model."""
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
X=sm.add_constant(der[['age']+BIN].values.astype(float)); y=der['amd_self'].values.astype(float)
m=sm.Logit(y,X).fit(disp=0)
names=['const','age']+BIN; conf=m.conf_int()
rows=[]
for i,nm in enumerate(names):
    rows.append(dict(term=nm,beta=round(float(m.params[i]),4),se=round(float(m.bse[i]),4),
        OR=round(float(np.exp(m.params[i])),4),OR_lo=round(float(np.exp(conf[i][0])),4),OR_hi=round(float(np.exp(conf[i][1])),4),
        p=round(float(m.pvalues[i]),5)))
apparent=roc_auc_score(y,m.predict(X))
# bootstrap optimism-corrected AUROC
rng=np.random.RandomState(SEED);opt=[]
for _ in range(500):
    idx=rng.randint(0,len(y),len(y)); yy=y[idx]
    if yy.sum() in (0,len(yy)):continue
    try:
        mb=sm.Logit(yy,X[idx]).fit(disp=0)
        opt.append(roc_auc_score(yy,mb.predict(X[idx]))-roc_auc_score(y,mb.predict(X)))
    except Exception: pass
corrected=float(apparent-np.mean(opt))
out=dict(derivation_n=int(len(y)),events=int(y.sum()),prev=round(float(y.mean()),4),
    apparent_auc=round(float(apparent),4),corrected_auc=round(corrected,4),
    brier=round(float(brier_score_loss(y,m.predict(X))),5),
    epv=round(float(y.sum()/len(BIN)),1),epv_by_df=round(float(y.sum()/(len(BIN)+1)),1),
    per_decade_OR=round(float(np.exp(m.params[names.index('age')]*10)),3),
    per_year_OR=round(float(np.exp(m.params[names.index('age')])),4),
    coefficients=rows)
json.dump(out,open(OUT/"sr_derivation.json","w"),indent=2)
print(json.dumps(out,indent=1)[:1800])
