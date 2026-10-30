"""Build NHIS (2016, 2017) harmonized external-validation dataset for Topic 1.
Analysis unit = person (sample adult). Outcome = self-reported AMD (VIM_MDEV==1).
Same 7 score predictors as NHANES; HTN/cholesterol/diabetes/obesity are ever-told / self-reported
here (vs measured in NHANES) -> intentional measurement-mode difference.
Codes verified against actual value distributions in _probe_nhis.py before use.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[v]="1"
import pandas as pd, numpy as np, json, tempfile

OUT = r"H:\眼科论文\题目4\analysis\output"
os.makedirs(OUT, exist_ok=True)
NH = r"H:\眼科论文\题目4\data\raw\nhis"
MISS = {7,9}

def load(y):
    df = pd.read_csv(os.path.join(NH, f"e{y}", "samadult.csv"))
    df.columns = [c.upper() for c in df.columns]
    for c in ['AGE_P','SEX','SMKSTAT2','BMI','HYPEV','CHLEV','DIBEV1','VIM_MDEV','WTFA_SA','PPSU','PSTRAT','REGION']:
        if c in df.columns: df[c]=pd.to_numeric(df[c],errors='coerce')
    o = pd.DataFrame()
    o['age']=df['AGE_P']
    o['male']=np.where(df['SEX'].notna(),(df['SEX']==1).astype(float),np.nan)
    sm=df['SMKSTAT2']
    o['smoke_current']=np.where(sm.notna(), sm.isin([1,2]).astype(float), np.nan)
    o['smoke_former']=np.where(sm.notna(), (sm==3).astype(float), np.nan)
    b=df['BMI']
    ob_valid=b.between(1000,6000)
    o['obesity']=np.where(ob_valid,(b>=3000).astype(float),np.nan)
    for src,var in [('HYPEV','hypertension'),('CHLEV','dyslipidemia')]:
        s=df[src].where(~df[src].isin(MISS))
        o[var]=np.where(s.notna(),(s==1).astype(float),np.nan)
    dia=df['DIBEV1'].where(~df['DIBEV1'].isin(MISS))  # 1=yes,2=no,3=borderline->no
    o['diabetes']=np.where(dia.notna(),(dia==1).astype(float),np.nan)
    amd=df['VIM_MDEV'].where(~df['VIM_MDEV'].isin(MISS))
    o['amd_self']=np.where(amd.notna(),(amd==1).astype(float),np.nan)
    o['weight']=df['WTFA_SA']; o['psu']=df['PPSU']; o['stratum']=df['PSTRAT']; o['region']=df['REGION']
    o['survey_year']=int(y)
    return o[(o['age']>=40)&o['amd_self'].notna()].copy()

frames=[load(y) for y in (2016,2017)]
d=pd.concat(frames,ignore_index=True)
preds=['age','male','smoke_former','smoke_current','obesity','hypertension','diabetes','dyslipidemia']
d['cc_pred']=d[preds].notna().all(axis=1)

tmp=tempfile.NamedTemporaryFile('w',delete=False,suffix='.csv',dir=OUT,newline='').name
d.to_csv(tmp,index=False); os.replace(tmp,os.path.join(OUT,"nhis_external.csv"))

summ={"n_40plus_valid_amd":int(len(d)),
 "by_year":{int(y):{"n":int((d.survey_year==y).sum()),
        "amd_prev":round(float(d.loc[d.survey_year==y,'amd_self'].mean()),5),
        "cc_n":int(((d.survey_year==y)&d.cc_pred).sum())} for y in (2016,2017)},
 "amd_prev_pooled":round(float(d['amd_self'].mean()),5),
 "exposure_prevalence_cc":{p:round(float(d.loc[d.cc_pred,p].mean()),4) for p in preds},
 "predictor_missing":{p:int(d[p].isna().sum()) for p in preds},
 "events_amd_cc":int(d.loc[d.cc_pred,'amd_self'].sum())}
with open(os.path.join(OUT,"nhis_build_summary.json"),"w") as f: json.dump(summ,f,indent=2)
print(json.dumps(summ,indent=2))
