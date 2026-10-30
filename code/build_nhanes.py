"""Build NHANES derivation analytic dataset for Topic 1 (frozen admission_record.md).

Lock BLAS/OMP threads BEFORE importing numpy (execution hygiene #7).
Analysis unit = person. Cycles 2005-06 (_D, SDDSRVYR=4) + 2007-08 (_E, =5).
Outcome primary Y = self-reported AMD (VIQ310==1). Objective ARM (OPXRET) added as
SEPARATE sensitivity columns (never pooled with self-report).
7 score predictors harmonized to be reproducible in NHIS: age, male, smoking(former/current),
obesity(BMI>=30), hypertension(measured >=140/90), diabetes(DIQ010==1), dyslipidemia(measured).
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[v] = "1"
import pyreadstat, pandas as pd, numpy as np, json, tempfile

RAW = r"H:\眼科论文\眼科2\data\raw\nhanes"
OUT = r"H:\眼科论文\题目4\analysis\output"
os.makedirs(OUT, exist_ok=True)

def read(f):
    df, _ = pyreadstat.read_xport(os.path.join(RAW, f))
    df.columns = [c.upper() for c in df.columns]
    return df

def cyc(suf):
    a, b = read(f"{suf}_D.xpt"), read(f"{suf}_E.xpt")
    return pd.concat([a, b], ignore_index=True)

demo = cyc("DEMO")[['SEQN','SDDSRVYR','RIDSTATR','RIAGENDR','RIDAGEYR','RIDRETH1','DMDEDUC2','INDFMPIR','WTMEC2YR','SDMVPSU','SDMVSTRA']]
viq  = cyc("VIQ")[['SEQN','VIQ071','VIQ090','VIQ310']]
diq  = cyc("DIQ")[['SEQN','DIQ010']].drop_duplicates('SEQN')
bpx  = cyc("BPX")[['SEQN','BPXSY2','BPXSY3','BPXSY4','BPXDI2','BPXDI3','BPXDI4']]
bmx  = cyc("BMX")[['SEQN','BMXBMI']].drop_duplicates('SEQN')
tc   = cyc("TCHOL")[['SEQN','LBXTC']].drop_duplicates('SEQN')
hdl  = cyc("HDL")[['SEQN','LBDHDD']].drop_duplicates('SEQN')
trg  = cyc("TRIGLY")[['SEQN','LBXTR','LBDLDL']].drop_duplicates('SEQN')
smq  = cyc("SMQ")[['SEQN','SMQ020','SMQ040']].drop_duplicates('SEQN')
# objective ARM worst-eye recodes (present in OPXRET_D/E)
opx  = cyc("OPXRET")[['SEQN','OPDUARM','OPDULARM']].drop_duplicates('SEQN')

d = demo.copy()
for r in [viq, diq, bpx, bmx, tc, hdl, trg, smq, opx]:
    d = d.merge(r, on='SEQN', how='left')

# numeric coercion
numcols = ['RIAGENDR','RIDAGEYR','RIDSTATR','SDDSRVYR','WTMEC2YR','SDMVPSU','SDMVSTRA','RIDRETH1','DMDEDUC2','INDFMPIR',
           'VIQ071','VIQ090','VIQ310','DIQ010','SMQ020','SMQ040','BMXBMI','LBXTC','LBDHDD','LBXTR','LBDLDL',
           'BPXSY2','BPXSY3','BPXSY4','BPXDI2','BPXDI3','BPXDI4','OPDUARM','OPDULARM']
for c in numcols:
    if c in d: d[c] = pd.to_numeric(d[c], errors='coerce')

# recode NHANES skip/don't-know codes to NaN for binary questionnaire items
for c in ['VIQ071','VIQ090','VIQ310','DIQ010','SMQ020','SMQ040']:
    d[c] = d[c].replace({7.0:np.nan, 9.0:np.nan, 77.0:np.nan, 99.0:np.nan})

# BP mean of readings 2-4
sbp = d[['BPXSY2','BPXSY3','BPXSY4']].mean(axis=1, skipna=True)
dbp = d[['BPXDI2','BPXDI3','BPXDI4']].mean(axis=1, skipna=True)
d['SBP_mean'] = sbp; d['DBP_mean'] = dbp

# --- predictors (frozen) ---
d['age'] = d['RIDAGEYR']
d['male'] = (d['RIAGENDR']==1).astype(float)
d.loc[d['RIAGENDR'].isna(),'male']=np.nan
# smoking: never ref / former / current
d['smoke_former'] = np.nan; d['smoke_current'] = np.nan
ever = d['SMQ020']==1
d.loc[d['SMQ020'].notna(),'smoke_former'] = (ever & (d['SMQ040']==3)).astype(float)
d.loc[d['SMQ020'].notna(),'smoke_current'] = (ever & d['SMQ040'].isin([1,2])).astype(float)
d['obesity'] = np.where(d['BMXBMI'].notna(), (d['BMXBMI']>=30).astype(float), np.nan)
d['hypertension'] = np.where(sbp.notna()|dbp.notna(), ((sbp>=140)|(dbp>=90)).astype(float), np.nan)
d['diabetes'] = np.where(d['DIQ010'].notna(), (d['DIQ010']==1).astype(float), np.nan)
sex = d['RIAGENDR']
hi_tc = d['LBXTC']>=240
lo_hdl = ((sex==1)&(d['LBDHDD']<40)) | ((sex==2)&(d['LBDHDD']<50))
hi_tg = d['LBXTR']>=200
hi_ldl = d['LBDLDL']>=160
dys = hi_tc | lo_hdl | hi_tg | hi_ldl
has_lipid = d[['LBXTC','LBDHDD','LBXTR','LBDLDL']].notna().any(axis=1)
d['dyslipidemia'] = np.where(has_lipid, dys.fillna(False).astype(float), np.nan)
# outcome + objective sensitivity
d['amd_self'] = np.where(d['VIQ310'].notna(), (d['VIQ310']==1).astype(float), np.nan)
d['cataract_surg'] = np.where(d['VIQ071'].notna(), (d['VIQ071']==1).astype(float), np.nan)
d['glaucoma_self'] = np.where(d['VIQ090'].notna(), (d['VIQ090']==1).astype(float), np.nan)
d['obj_anyARM']  = np.where(d['OPDUARM'].notna(),  (d['OPDUARM']>=1).astype(float), np.nan)
d['obj_lateARM'] = np.where(d['OPDULARM'].notna(), (d['OPDULARM']==1).astype(float), np.nan)

# eligibility: age>=40 with valid self-reported AMD outcome
d = d[(d['age']>=40) & d['amd_self'].notna()].copy()
preds = ['age','male','smoke_former','smoke_current','obesity','hypertension','diabetes','dyslipidemia']
d['cc_pred'] = d[preds].notna().all(axis=1)

outcols = ['SEQN','SDDSRVYR','RIDSTATR','WTMEC2YR','SDMVPSU','SDMVSTRA','RIDRETH1','DMDEDUC2','INDFMPIR',
           'age','male','smoke_former','smoke_current','obesity','hypertension','diabetes','dyslipidemia',
           'amd_self','cataract_surg','glaucoma_self','obj_anyARM','obj_lateARM',
           'SBP_mean','DBP_mean','BMXBMI','LBXTC','LBDHDD','LBXTR','LBDLDL','cc_pred']
d = d[outcols]

# atomic write
tmp = tempfile.NamedTemporaryFile('w', delete=False, suffix='.csv', dir=OUT, newline='').name
d.to_csv(tmp, index=False)
os.replace(tmp, os.path.join(OUT, "nhanes_derivation.csv"))

summ = {
  "n_age40_valid_amd": int(len(d)),
  "cycle_counts": d['SDDSRVYR'].value_counts().to_dict(),
  "n_complete_case": int(d['cc_pred'].sum()),
  "prevalence_amd_self": round(float(d['amd_self'].mean()),5),
  "prevalence_amd_self_cc": round(float(d.loc[d['cc_pred'],'amd_self'].mean()),5),
  "prevalence_obj_anyARM": round(float(d['obj_anyARM'].mean()),5),
  "prevalence_obj_lateARM": round(float(d['obj_lateARM'].mean()),5),
  "n_valid_obj_anyARM": int(d['obj_anyARM'].notna().sum()),
  "events_amd_self_cc": int(d.loc[d['cc_pred'],'amd_self'].sum()),
  "predictor_missing": {p: int(d[p].isna().sum()) for p in preds},
  "exposure_prevalence_cc": {p: round(float(d.loc[d['cc_pred'],p].mean()),4) for p in preds},
}
with open(os.path.join(OUT,"nhanes_build_summary.json"),"w") as f:
    json.dump(summ,f,indent=2)
print(json.dumps(summ,indent=2))
print("n rows written:", len(d))
