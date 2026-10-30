"""Self-contained published risk equation (Harrell 3-df RCS on age, knots [41,53,65,80]).
No patsy in the shipped calculator -> fully reproducible from the coefficient list alone.
Writes analysis/output/model_equation.json and manuscript/risk_calculator.py, and validates.
"""
import os, numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
OUT=r"H:\眼科论文\题目4\analysis\output"
KNOTS=[41.0,53.0,65.0,80.0]
BIN=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
def rcs(x):
    t1,t2,t3,t4=KNOTS; x=np.asarray(x,float); d=t4-t3
    cp=lambda t: np.maximum(x-t,0)**3
    r1=cp(t1)-cp(t3)*(t4-t1)/d+cp(t4)*(t3-t1)/d
    r2=cp(t2)-cp(t3)*(t4-t2)/d+cp(t4)*(t3-t2)/d
    return x,r1,r2
der=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv")); der=der[der.cc_pred==True].reset_index(drop=True)
a,r1,r2=rcs(der['age'].values)
D=pd.DataFrame(dict(age=a,age_r1=r1,age_r2=r2,**{b:der[b].values for b in BIN}))
Xd=sm.add_constant(D); yd=der['amd_self'].values.astype(float)
m=sm.Logit(yd,Xd.values).fit(disp=0)
names=['const']+list(D.columns); coef=dict(zip(names,[round(float(v),5) for v in m.params]))
auc=round(float(roc_auc_score(yd,m.predict(Xd.values))),4)
lp_ref="logit(p) = "+ " + ".join([f"({coef['const']})"]+[f"({coef['age']})*age",f"({coef['age_r1']})*R1(age)",f"({coef['age_r2']})*R2(age)"]+[f"({coef[b]})*{b}" for b in BIN])

def prob(age,male=0,smoke_former=0,smoke_current=0,obesity=0,hypertension=0,diabetes=0,dyslipidemia=0):
    A,R1,R2=rcs([age]); lp=coef['const']+coef['age']*A[0]+coef['age_r1']*R1[0]+coef['age_r2']*R2[0]
    for name,val in zip(BIN,[male,smoke_former,smoke_current,obesity,hypertension,diabetes,dyslipidemia]):
        lp+=coef[name]*val
    return round(float(lp),4), round(float(1/(1+np.exp(-lp))),4)
hi=prob(70,male=1,smoke_former=1,obesity=1,hypertension=1,dyslipidemia=1)
lo=prob(55)
agecurve=[[a,prob(a)[1]] for a in range(45,86,5)]
eq=dict(model="Logistic, outcome = self-reported AMD, NHANES 2005-08 complete cases (n=%d, %d events); derivation apparent AUROC %.4f"%(len(der),int(yd.sum()),auc),
    age_basis="Harrell restricted cubic spline, 3 df, knots=[41,53,65,80] (R1,R2 per Harrell truncated-cubic definition; see risk_calculator.py)",
    reference_covariates="male=0 (female), non-smoker, non-obese, normotensive, non-diabetic, normal lipids",
    coefficients=coef, lp_formula=lp_ref,
    worked_examples={"high_risk_70yM_formerSmoker_obese_HTN_highchol":{"lp":hi[0],"p":hi[1]},
                     "low_risk_55yF_healthy":{"lp":lo[0],"p":lo[1]}},
    age_only_reference_risk_curve=agecurve)
json.dump(eq,open(os.path.join(OUT,"model_equation.json"),"w"),indent=2)

# ---- self-contained calculator (no patsy) ----
calc='''"""Frozen AMD clinical risk score (Topic 1). Pure-python, no patsy required.
logit(p)=c0 + b_age*age + b_r1*R1(age) + b_r2*R2(age) + sum b_i * x_i
Knots (Harrell RCS, 3 df) = %s
Coefficients: %s
"""
import math
KNOTS=%s
C=%s
def _cp(x,t): return max(x-t,0.0)**3
def _rcs(x):
    t1,t2,t3,t4=KNOTS; d=t4-t3
    r1=_cp(x,t1)-_cp(x,t3)*(t4-t1)/d+_cp(x,t4)*(t3-t1)/d
    r2=_cp(x,t2)-_cp(x,t3)*(t4-t2)/d+_cp(x,t4)*(t3-t2)/d
    return x,r1,r2
def amd_risk_logit(age,male=0,smoke_former=0,smoke_current=0,obesity=0,hypertension=0,diabetes=0,dyslipidemia=0):
    A,R1,R2=_rcs(age)
    return C['const']+C['age']*A+C['age_r1']*R1+C['age_r2']*R2+C['male']*male+C['smoke_former']*smoke_former+C['smoke_current']*smoke_current+C['obesity']*obesity+C['hypertension']*hypertension+C['diabetes']*diabetes+C['dyslipidemia']*dyslipidemia
def amd_risk_prob(**kw):
    lp=amd_risk_logit(**kw); return 1.0/(1.0+math.exp(-lp))
if __name__=="__main__":
    print("70y M (former smoker, obese, HTN, high-chol, non-diab):",amd_risk_prob(age=70,male=1,smoke_former=1,obesity=1,hypertension=1,dyslipidemia=1))
    print("55y F healthy:",amd_risk_prob(age=55))
'''%(KNOTS,coef,KNOTS,coef)
open(r"H:\眼科论文\题目4\manuscript\risk_calculator.py","w").write(calc)
print("apparent AUROC:",auc); print("high70M:",hi,"low55F:",lo)
print("LP formula:",lp_ref[:120])
print("calc self-check:")
