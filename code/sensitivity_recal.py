"""Sensitivity: recalibration tiers + design-weighted replication for Topic 1.
Quantifies how much of the observed calibration drift is fixed by intercept-only vs
intercept+slope recalibration, and whether weighted vs unweighted derivation changes conclusions.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[v]="1"
import numpy as np, pandas as pd, json, patsy, statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit
import statsmodels.api as sm_api
from statsmodels.genmod.generalized_linear_model import GLM
from sklearn.metrics import roc_auc_score, brier_score_loss
SEED=20261028
OUT=r"H:\眼科论文\题目4\analysis\output"
PRED=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
RHS="cr(age, df=3) + "+" + ".join(PRED)
der=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv")); der=der[der.cc_pred==True].copy()
ext=pd.read_csv(os.path.join(OUT,"nhis_external.csv")); ext=ext[ext.cc_pred==True].copy()
yder,Xder=patsy.dmatrices("amd_self ~ "+RHS,der,return_type="dataframe")
di=Xder.design_info
m=Logit(yder.values.ravel(),Xder.values).fit(disp=0); coef=m.params
# frozen LP on NHIS pooled
Xext=patsy.build_design_matrices([di],ext,return_type="dataframe")[0]
lp=Xext.values@coef; y=ext['amd_self'].values.astype(float); prev=y.mean()
def metrics(y,p):
    return dict(auc=round(float(roc_auc_score(y,p)),4),brier=round(float(brier_score_loss(y,p)),5),
                mean_pred=round(float(p.mean()),4))
res={}
# as-is
p_asis=1/(1+np.exp(-lp)); res['asis']=metrics(y,p_asis)
# intercept-only recalibration: shift LP by calib_intercept (slope stays 1)
Xi=sm.add_constant(np.ones(len(y))*0)  # placeholder
fint=Logit(y,np.ones((len(y),1)),offset=lp).fit(disp=0); a=fint.params[0]
p_int=1/(1+np.exp(-(lp+a))); res['intercept_only']=metrics(y,p_int); res['intercept_only']['shift']=round(float(a),4)
# intercept+slope recalibration
Xs=sm.add_constant(lp); fs=Logit(y,Xs).fit(disp=0); b0,b1=fs.params
p_full=1/(1+np.exp(-(b0+b1*lp))); res['intercept_slope']=metrics(y,p_full)
res['recalibration_coefs']=dict(intercept_only_a=round(float(a),4),full_b0=round(float(b0),4),full_b1=round(float(b1),4))
# design-weighted derivation model (GLM Binomial with weight=WTMEC2YR/2 combined-> use as freq approx)
w=der['WTMEC2YR'].values/2.0
try:
    gw=GLM(yder.values.ravel(),Xder.values,family=sm_api.families.Binomial(),freq_weights=w).fit()
    lpw=patsy.build_design_matrices([di],ext,return_type="dataframe")[0].values@gw.params
    res['weighted_derivation']=dict(auc_nhis=round(float(roc_auc_score(y,1/(1+np.exp(-lpw)))),4),
                                    note="freq_weights approx of design weighting")
except Exception as e:
    res['weighted_derivation']=f"err:{e}"
json.dump(res,open(os.path.join(OUT,"sensitivity_recal.json"),"w"),indent=2)
print(json.dumps(res,indent=2))
