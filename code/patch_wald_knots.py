"""Patch: joint Wald chi2 for age spline terms + extract RCS knots; update harden_results.json."""
import os, numpy as np, pandas as pd, patsy, json
from scipy import stats
from statsmodels.discrete.discrete_model import Logit
OUT=r"H:\眼科论文\题目4\analysis\output"
BIN=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
der=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv")); der=der[der.cc_pred==True].reset_index(drop=True)
y,X=patsy.dmatrices("amd_self ~ cr(age, df=3) "+" ".join(["+ "+b for b in BIN]),der,return_type="dataframe")
m=Logit(y.values.ravel(),X.values).fit(disp=0)
cols=list(X.design_info.column_names)
age_idx=[cols.index(c) for c in cols if c.startswith("cr(age")]
b=np.asarray(m.params)[age_idx]; V=np.asarray(m.cov_params())[np.ix_(age_idx,age_idx)]
W=float(b@np.linalg.solve(V,b)); df=len(age_idx); p=float(stats.chi2.sf(W,df))
# knots extraction from design_info
knots=None
try:
    for ti in X.design_info.term_transforms_info:
        st=ti.state
        ks=getattr(st,'knots',None)
        if ks is not None:
            knots=[round(float(k),2) for k in np.ravel(ks)]; break
except Exception as e:
    print("knots err",e)
# fallback: replicate patsy default knots for df=3 (quantiles)
if not knots:
    q=[0.05,0.5,0.95] if False else None
    # patsy cr default: knots at quantiles (0,1)/(n+1)... use numpy to match len
    knots=None
print("AGE joint Wald chi2=%.2f df=%d p=%.3g"%(W,df,p))
print("knots:",knots, "| design_info repr:")
try:
    print(X.design_info.term_transforms_info[0].state)
except Exception as e:
    print(e)
h=json.load(open(os.path.join(OUT,"harden_results.json")))
h['age_wald']=dict(statistic=round(W,2),df=df,p=p,method="Wald chi-square (cr(age) terms jointly = 0)")
if knots: h['knots']=knots
json.dump(h,open(os.path.join(OUT,"harden_results.json"),"w"),indent=2)
print("patched age_wald (and knots if found)")
