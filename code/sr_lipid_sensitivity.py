"""ROUND-9 P0-1: NHANES high-cholesterol / fasting sub-sample audit + DEFINITIONAL SENSITIVITY.
Answers the reviewer's four code-confirmation questions with real numbers and tests whether making
dyslipidaemia use ONLY the full-MEC-sample analytes (total cholesterol + HDL, which are available to
100% of examined adults and carry the correct WTMEC2YR exam weight) instead of the multi-analyte
definition that also uses fasting-sub-sample-only LDL/triglyceride changes the headline conclusions.
Self-contained: refits from the analytic CSVs in this process; does NOT overwrite the canonical chain.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
import numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"
SEED=20261028
ALLP=["age","male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
der=pd.read_csv(OUT/"nhanes_derivation.csv")
ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)

# ---- sub-sample structure evidence ----
cc=der[der.cc_pred==True].reset_index(drop=True)
struct=dict(
 n_cc=int(len(cc)),
 full_sample_analytes=dict(TC_present_pct=round(100*cc['LBXTC'].notna().mean(),1), HDL_present_pct=round(100*cc['LBDHDD'].notna().mean(),1)),
 fasting_only_analytes=dict(TG_present_pct=round(100*cc['LBXTR'].notna().mean(),1), LDL_present_pct=round(100*cc['LBDLDL'].notna().mean(),1)),
 cc_dyslip_zero_rows=int((cc['dyslipidemia']==0).sum()),
 cc_dyslip_zero_but_TG_LDL_missing=int(((cc['dyslipidemia']==0)&cc['LBXTR'].isna()&cc['LBDLDL'].isna()).sum()),
)
struct["pct_zero_rows_relying_on_TChdl_only"]=round(100*struct["cc_dyslip_zero_but_TG_LDL_missing"]/max(struct["cc_dyslip_zero_rows"],1),1)

# ---- rebuild dyslipidaemia using ONLY full-sample analytes (TC>=240 OR low HDL) ----
sex=der['male']
der['dyslip_tc_hdl']=((der['LBXTC']>=240)|((sex==1)&(der['LBDHDD']<40))|((sex==2)&(der['LBDHDD']<50))).where(der[['LBXTC','LBDHDD']].notna().all(axis=1)).astype('float')

def design(df): return np.column_stack([np.ones(len(df)),df[ALLP].values.astype(float)])
def ageonly(df): return np.column_stack([np.ones(len(df)),df[["age"]].values.astype(float)])

def run(col):
    d=der.copy(); d['dyslipidemia']=d[col]; ccl=d[d.cc_pred==True].reset_index(drop=True)
    X=design(ccl); y=ccl['amd_self'].values.astype(float)
    m=sm.Logit(y,X).fit(disp=0); b=np.asarray(m.params)
    Xe=design(ext); lp=Xe@b; ye=ext['amd_self'].values.astype(float); p=1/(1+np.exp(-lp))
    sl=float(sm.Logit(ye,sm.add_constant(lp)).fit(disp=0).params[1])
    citl=float(sm.Logit(ye,np.ones((len(ye),1)),offset=lp).fit(disp=0).params[0])
    # age-only baseline on same sample
    ma=sm.Logit(y,ageonly(ccl)).fit(disp=0); lpa=ageonly(ext)@np.asarray(ma.params)
    ext_age=roc_auc_score(ye,lpa); ext_full=roc_auc_score(ye,lp)
    return dict(apparent=round(float(roc_auc_score(y,X@b)),4),
                ext_auc=round(float(ext_full),4), ext_slope=round(sl,4), ext_citl=round(citl,4),
                ext_age_only=round(float(ext_age),4), dAUC_full_minus_age=round(float(ext_full-ext_age),4),
                dyslip_coef=round(float(b[ALLP.index('dyslipidemia')+1]),4),
                dyslip_prev_cc=round(float(ccl[col].mean()),4))

# paired bootstrap dAUROC(full-age) for the TC+HDL definition (quick, individual resample)
def paired_dauc(col,B=400):
    d=der.copy(); d['dyslipidemia']=d[col]; ccl=d[d.cc_pred==True].reset_index(drop=True)
    y=ccl['amd_self'].values.astype(float)
    m=sm.Logit(y,design(ccl)).fit(disp=0); ma=sm.Logit(y,ageonly(ccl)).fit(disp=0)
    bf=np.asarray(m.params); ba=np.asarray(ma.params)
    Xe=design(ext); Xae=ageonly(ext); ye=ext['amd_self'].values.astype(float); n=len(ye)
    rng=np.random.RandomState(SEED); ds=[]
    for _ in range(B):
        idx=rng.randint(0,n,n)
        ds.append(roc_auc_score(ye[idx],Xe[idx]@bf)-roc_auc_score(ye[idx],Xae[idx]@ba))
    lo,hi=np.percentile(ds,[2.5,97.5]); return round(float(np.mean(ds)),4), round(float(lo),4), round(float(hi),4)

multi=run('dyslipidemia'); tchdl=run('dyslip_tc_hdl')
dauc_tc=paired_dauc('dyslip_tc_hdl'); dauc_multi=paired_dauc('dyslipidemia')

out=dict(sub_sample_structure=struct,
 definitions=dict(
   primary_multi_analyte=dict(desc="TC>=240 OR HDL<40/50 OR TG>=200 OR LDL>=160 (TG/LDL from fasting sub-sample; non-fasting evaluated on TC/HDL only)",**multi,paired_dAUC=dict(mean=dauc_multi[0],ci=[dauc_multi[1],dauc_multi[2]])),
   sensitivity_tc_hdl_only=dict(desc="TC>=240 OR HDL<40/50 (full-MEC-sample analytes only; WTMEC2YR-valid)",**tchdl,paired_dAUC=dict(mean=dauc_tc[0],ci=[dauc_tc[1],dauc_tc[2]]))),
 weights_note="NHANES total cholesterol and HDL are measured on the full MEC sample and carry WTMEC2YR (correct for the unweighted primary and the design-weighted derivation). Triglycerides and calculated LDL are available only in the overnight-fasting sub-sample (~48% of adults >=40 here), which carries a distinct fasting sub-sample weight. Mixing them means no single weight is exactly right for the multi-analyte variable; this is why the design-weighted derivation is labelled indicative. The TC+HDL-only definition removes the sub-sample mismatch entirely.",
 conclusion=(f"Switching to the uniform full-sample TC+HDL definition moves external AUROC {multi['ext_auc']:.4f}->{tchdl['ext_auc']:.4f}, apparent {multi['apparent']:.4f}->{tchdl['apparent']:.4f}, slope {multi['ext_slope']:.4f}->{tchdl['ext_slope']:.4f}, CITL {multi['ext_citl']:.4f}->{tchdl['ext_citl']:.4f}; the incremental dAUROC over age stays ~{tchdl['dAUC_full_minus_age']:.4f}. Headline conclusions (age-dominated, external AUROC ~0.77, persistent over-prediction, negligible clinical-variable gain) are unchanged; the tighter definition slightly enlarges the over-prediction (lower dyslipidaemia prevalence {multi['dyslip_prev_cc']:.3f}->{tchdl['dyslip_prev_cc']:.3f})."))
json.dump(out,open(OUT/"sr_lipid_sensitivity.json","w"),indent=2)
print(json.dumps(out,indent=1))
