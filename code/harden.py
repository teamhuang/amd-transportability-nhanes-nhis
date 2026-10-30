"""HARDENING re-analysis for Topic 1 (addresses reviewer round). Deterministic, single source.
Re-derives with EXPLICIT spline knots (transparent, transportable), and adds:
 nested model comparison (age-only/+sex/full), DeltaAUC CI, common age-range (40-69),
 measurement-mode-clean subset, design-weighted estimation + PSU-cluster bootstrap CIs,
 temporal recalibration (2016<->2017), missingness, subgroups, knot sensitivity,
 objective-endpoint definitions + gradeability counts + age-only baseline,
 full risk equation + knots + worked example + EPV-by-df, corrected Holm(7) + joint age Wald,
 prespecified DCA thresholds + NNE.
All outputs -> analysis/output/harden_*.json and model_equation.json.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[v]="1"
import numpy as np, pandas as pd, json, patsy, tempfile
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Binomial
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_auc_score, brier_score_loss
SEED=20261028
OUT=r"H:\眼科论文\题目4\analysis\output"
def J(o,n):
    tmp=tempfile.NamedTemporaryFile('w',delete=False,suffix='.json',dir=OUT,newline='').name
    json.dump(o,open(tmp,'w'),indent=2); os.replace(tmp,os.path.join(OUT,n))
BIN=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
CLEAN=["male","smoke_former","smoke_current","diabetes"]           # drop mode-divergent HTN/lipid & measured-vs-self obesity
MODE_DIV=["obesity","hypertension","dyslipidemia"]

der=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv")); der=der[der.cc_pred==True].reset_index(drop=True)
ext=pd.read_csv(os.path.join(OUT,"nhis_external.csv")); ext=ext[ext.cc_pred==True].reset_index(drop=True)
from rcs_basis import add_basis
der=add_basis(der); ext=add_basis(ext)

# ---- spline: standard patsy cr(df) (df=3 => 3 columns). Extract actual knots for transparency ----
def get_knots(Xdf):
    try:
        for t in Xdf.design_info.term_transforms_info:
            st=t.state
            try:
                return [round(float(k),2) for k in st.knots]
            except Exception:
                pass
        # CRSplineState may expose knots as attribute list
        return None
    except Exception:
        return None
def rhs(binp,knotdf=3):
    spline = "age_lin + age_r1 + age_r2" if knotdf>=3 else ("age_lin" if knotdf==1 else "age_lin + age_r1 + age_r2")
    return spline + ("".join([" + "+b for b in binp]) if binp else "")

def fit(rhsstr,data):
    y,X=patsy.dmatrices("amd_self ~ "+rhsstr,data,return_type="dataframe")
    m=Logit(y.values.ravel(),X.values).fit(disp=0)
    return m,X,y

def lp_apply(model_info,newdata):
    di=model_info['design_info']; coef=model_info['coef']
    Xn=patsy.build_design_matrices([di],newdata,return_type="dataframe")[0]
    return Xn.values@coef, np.array(Xn.columns.tolist())

# ---------- nested models (full transportable pipeline) ----------
def transport_metrics(y,lp,weight=None):
    p=1/(1+np.exp(-lp))
    auc=roc_auc_score(y,lp)
    slope=float(sm.Logit(y.astype(float),sm.add_constant(lp)).fit(disp=0).params[1])
    ci0=float(sm.Logit(y.astype(float),np.ones((len(y),1)),offset=lp).fit(disp=0).params[0])
    return dict(auc=round(float(auc),4),brier=round(float(brier_score_loss(y,p)),5),
                calib_slope=round(slope,4),calib_intercept=round(ci0,4),
                mean_pred=round(float(p.mean()),4),obs=round(float(np.mean(y)),4))

models={}
specs={"age_only":rhs([],3),"age_sex":rhs(["male"],3),"full":rhs(BIN,3),"mode_clean":rhs(CLEAN,3)}
frozen={}
for name,r in specs.items():
    m,X,y=fit(r,der)
    frozen[name]=dict(design_info=X.design_info,coef=m.params)
    dm={'converged':bool(m.mle_retvals['converged'])}
    ext_m=transport_metrics(ext['amd_self'].values.astype(float),*lp_apply(frozen[name],ext) )
    models[name]=dict(derivation=dict(n=int(len(y)),events=int(y.values.sum()),
                    apparent_auc=round(float(roc_auc_score(y.values,m.predict(X.values))),4),
                    terms=X.shape[1], **dm),
                      external_pooled=ext_m)
# DeltaAUC full vs age_only on external (bootstrap, individual)
lp_full,_=lp_apply(frozen['full'],ext); lp_age,_=lp_apply(frozen['age_only'],ext)
yext=ext['amd_self'].values.astype(float)
rng=np.random.RandomState(SEED); deltas=[]
for _ in range(1000):
    idx=rng.randint(0,len(yext),len(yext)); yy=yext[idx]
    if yy.sum()==0 or yy.sum()==len(yy): continue
    deltas.append(roc_auc_score(yy,lp_full[idx])-roc_auc_score(yy,lp_age[idx]))
dlo,dhi=np.percentile(deltas,[2.5,97.5])
models['delta_auc_full_minus_ageonly']=dict(value=round(float(roc_auc_score(yext,lp_full)-roc_auc_score(yext,lp_age)),4),
    ci=[round(float(dlo),4),round(float(dhi),4)],nboot=len(deltas))

# ---------- common age range 40-69 ----------
e69=ext[(ext.age>=40)&(ext.age<=69)]; 
lp69,_=lp_apply(frozen['full'],e69)
models['age40_69']=transport_metrics(e69['amd_self'].values.astype(float),lp69)
d69=der[(der.age>=40)&(der.age<=69)]

# ---------- knot / functional-form sensitivity (linear vs 3-df RCS) on external AUC ----------
knot_sens={}
for kd,kdname in [(1,"linear_age"),(3,"rcs_3df")]:
    m,X,y=fit(rhs(BIN,kd),der)
    fi=dict(design_info=X.design_info,coef=m.params)
    lpf,_=lp_apply(fi,ext)
    knot_sens[kdname]=round(float(roc_auc_score(yext,lpf)),4)
models['knot_sensitivity']=knot_sens

# ---------- design-weighted estimation + PSU cluster bootstrap ----------
w_der=der['WTMEC2YR'].values/2.0   # MEC 2-year weights pooled over 2 cycles -> /2
yw,_=patsy.dmatrices("amd_self ~ "+rhs(BIN,3),der,return_type="dataframe")
gw=GLM(yw.values.ravel(),_ if (_ := None) else None,family=Binomial()) if False else None
# weighted prevalence
def wprev(v,w): 
    v=np.asarray(v,float); w=np.asarray(w,float); ok=~np.isnan(v)
    return round(float(np.nansum(v[ok]*w[ok])/np.nansum(w[ok])),4)
wprevs=dict(nhanes_amd_w=wprev(der['amd_self'],der['WTMEC2YR']),
            nhis_amd_w=wprev(ext['amd_self'],ext['weight']/2.0),
            nhis_htn_w=wprev(ext['hypertension'],ext['weight']/2.0),
            nhanes_htn_w=wprev(der['hypertension'],der['WTMEC2YR']))
try:
    gw=GLM(yw.values.ravel(),X if False else patsy.dmatrix(rhs(BIN,3),der,return_type="dataframe").values,
           family=Binomial(),freq_weights=w_der).fit()
    wcoef=gw.params
    weighted_note="weighted point estimates via GLM Binomial freq_weights"
except Exception as e:
    wcoef=None; weighted_note=f"weighted GLM failed: {e}"
# PSU-cluster bootstrap CI for external AUROC/slope/intercept (unweighted frozen model)
extc=ext.copy(); extc['psu_key']=extc['survey_year'].astype(str)+"_"+extc['stratum'].astype(str)+"_"+extc['psu'].astype(str)
groups=extc.groupby('psu_key')
gids=list(groups.groups.keys())
rng=np.random.RandomState(SEED); boot=[]
for _ in range(200):
    pick=rng.choice(len(gids),len(gids),replace=True)
    sub=pd.concat([groups.get_group(gids[i]) for i in pick])
    lps,_=lp_apply(frozen['full'],sub); ys=sub['amd_self'].values.astype(float)
    if ys.sum()==0 or ys.sum()==len(ys): continue
    try:
        sl=float(sm.Logit(ys,sm.add_constant(lps)).fit(disp=0).params[1])
        it=float(sm.Logit(ys,np.ones((len(ys),1)),offset=lps).fit(disp=0).params[0])
        boot.append((roc_auc_score(ys,lps),sl,it))
    except Exception: pass
b=np.array(boot)
psu_ci=dict(n_rep=len(boot),
    auc_ci=[round(float(np.percentile(b[:,0],2.5)),4),round(float(np.percentile(b[:,0],97.5)),4)],
    slope_ci=[round(float(np.percentile(b[:,1],2.5)),4),round(float(np.percentile(b[:,1],97.5)),4)],
    intercept_ci=[round(float(np.percentile(b[:,2],2.5)),4),round(float(np.percentile(b[:,2],97.5)),4)])
design=dict(weighted_prevalences=wprevs,weighted_note=weighted_note,
           weighted_betas=(None if wcoef is None else {c:round(float(wcoef[i]),4) for i,c in enumerate(patsy.dmatrix(rhs(BIN,3),der,return_type="dataframe").design_info.column_names)}),
           psu_cluster_bootstrap=psu_ci)

# ---------- temporal recalibration 2016<->2017 ----------
def recalib_intercept(source_y,source_lp):
    return float(sm.Logit(source_y.astype(float),np.ones((len(source_y),1)),offset=source_lp).fit(disp=0).params[0])
e16=ext[ext.survey_year==2016]; e17=ext[ext.survey_year==2017]
lp16,_=lp_apply(frozen['full'],e16); lp17,_=lp_apply(frozen['full'],e17)
a16=recalib_intercept(e16['amd_self'].values,lp16); a17=recalib_intercept(e17['amd_self'].values,lp17)
def calib_with(y,lp,a):
    p=1/(1+np.exp(-(lp+a)))
    return dict(brier=round(float(brier_score_loss(y,p)),5),mean_pred=round(float(p.mean()),4))
temporal=dict(
  a_from2016=round(a16,4),a_from2017=round(a17,4),
  naive_2017=calib_with(e17['amd_self'].values.astype(float),lp17,0.0),
  recal2016_on2017=calib_with(e17['amd_self'].values.astype(float),lp17,a16),
  naive_2016=calib_with(e16['amd_self'].values.astype(float),lp16,0.0),
  recal2017_on2016=calib_with(e16['amd_self'].values.astype(float),lp16,a17),
  auroc_frozen_2016=round(float(roc_auc_score(e16['amd_self'],lp16)),4),
  auroc_frozen_2017=round(float(roc_auc_score(e17['amd_self'],lp17)),4))

# ---------- missingness ----------
der_all=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv"))
miss_nhanes={p:round(float(der_all[p].isna().mean()),4) for p in ["obesity","hypertension","dyslipidemia","diabetes","smoke_former","smoke_current"]}
miss_nhanes["any_pred_missing_pct"]=round(float((~der_all.cc_pred.astype(bool)).mean()),4)
miss_nhanes["included_amd_prev"]=round(float(der[der.cc_pred==True]['amd_self'].mean()),4)
miss_nhanes["all_age40_amd_prev"]=round(float(der_all['amd_self'].mean()),4)
missingness=dict(nhanes=miss_nhanes)
ext_all=pd.read_csv(os.path.join(OUT,"nhis_external.csv"))
missingness["nhis_any_pred_missing_pct"]=round(float((~ext_all.cc_pred.astype(bool)).mean()),4)

# ---------- subgroups (NHIS, full model) ----------
lp_ext_all,_=lp_apply(frozen['full'],ext)
bands=[("40-54",40,54),("55-64",55,64),("65-74",65,74),("75+",75,200)]
subgroups=[]
for nm,lo,hi in bands:
    msk=(ext.age>=lo)&(ext.age<=hi)
    if msk.sum()>200:
        t=transport_metrics(ext['amd_self'].values[msk].astype(float),lp_ext_all[msk]); t['subgroup']=f"age {nm}"; subgroups.append(t)
for sx,lab in [(1,"male"),(0,"female")]:
    msk=ext['male'].values==sx
    if msk.sum()>200:
        t=transport_metrics(ext['amd_self'].values[msk].astype(float),lp_ext_all[msk]); t['subgroup']=f"sex {lab}"; subgroups.append(t)

# ---------- full coefficient table + Holm(7) + joint age Wald ----------
m_full,X_full,y_full=fit(rhs(BIN,3),der)
from rcs_basis import KNOTS as kn
cols=list(X_full.design_info.column_names)
pdict=dict(zip(cols,m_full.pvalues))
p7=[pdict[b] for b in BIN]
rej,holm,_,_=multipletests(np.array(p7),method='holm')
coef_table=[]
for i,c in enumerate(cols):
    conf=m_full.conf_int()[i]
    coef_table.append(dict(term=c,beta=round(float(m_full.params[i]),5),se=round(float(m_full.bse[i]),5),
        OR=round(float(np.exp(m_full.params[i])),4),OR_lo=round(float(np.exp(conf[0])),4),OR_hi=round(float(np.exp(conf[1])),4),
        p_raw=round(float(m_full.pvalues[i]),5)))
for j,b in enumerate(BIN): coef_table[[t['term'] for t in coef_table].index(b)]['p_holm7']=round(float(holm[j]),5)
# joint Wald chi2 for the age basis (age_lin, age_r1, age_r2)
from scipy import stats as _st
agecols=["age_lin","age_r1","age_r2"]
ai=[cols.index(c) for c in agecols if c in cols]
bav=np.asarray(m_full.params)[ai]; Vav=np.asarray(m_full.cov_params())[np.ix_(ai,ai)]
Wald=float(bav@np.linalg.solve(Vav,bav)); dfw=len(ai); pwald=float(_st.chi2.sf(Wald,dfw))
age_wald=dict(statistic=round(Wald,2),df=dfw,p=pwald,method="joint Wald chi-square (age basis terms jointly = 0)")
EPV=dict(events=int(y_full.values.sum()),predictor_terms=len(cols)-1,
         EPV_by_predictors=round(float(y_full.values.sum()/len(BIN)),1),
         EPV_by_df=round(float(y_full.values.sum()/(len(cols)-1)),1))

# ---------- worked example ----------
ex=pd.DataFrame({b:[0.0] for b in BIN})
ex['age']=70.0; ex['male']=1.0; ex['smoke_former']=1.0; ex['smoke_current']=0.0
ex['obesity']=1.0; ex['hypertension']=1.0; ex['diabetes']=0.0; ex['dyslipidemia']=1.0
ex=add_basis(ex)
lpe,_=lp_apply(dict(design_info=X_full.design_info,coef=m_full.params),ex)
worked=dict(desc="70y male, former smoker, obese, hypertensive, non-diabetic, high-cholesterol",
            linear_predictor=round(float(lpe[0]),4),probability=round(float(1/(1+np.exp(-lpe[0]))),4))

# ---------- objective endpoint: definitions, gradeability, age-only baseline, calibration ----------
der_raw=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv"))
valid_arm=der_raw['obj_anyARM'].notna().sum()
obj=dict(n_valid_worseye_grading=int(valid_arm),
         n_age40_cc=int(len(der)),
         pct_ungradeable_or_missing=round(float((len(der)- (der['obj_anyARM'].notna().sum()) )/len(der)),4))
obj_metrics={}
for ey,lab in [('obj_anyARM','any-ARM'),('obj_lateARM','late-AMD')]:
    s=der[der[ey].notna()]
    yv=s[ey].values.astype(float)
    if yv.sum()>=10:
        lf,_=lp_apply(frozen['full'],s); la,_=lp_apply(frozen['age_only'],s)
        obj_metrics[lab]=dict(n=int(len(s)),events=int(yv.sum()),prev=round(float(yv.mean()),5),
            auc_full=round(float(roc_auc_score(yv,lf)),4),auc_ageonly=round(float(roc_auc_score(yv,la)),4),
            calib_full=transport_metrics(yv,lf))
# concordance (self vs any-ARM) with corrected labels
comp=der[[ 'amd_self','obj_anyARM']].dropna()
conc=dict(n=int(len(comp)),fundus_pos_self_pos=int(((comp['amd_self']==1)&(comp['obj_anyARM']==1)).sum()),
          fundus_pos_self_neg=int(((comp['amd_self']==0)&(comp['obj_anyARM']==1)).sum()),
          fundus_neg_self_pos=int(((comp['amd_self']==1)&(comp['obj_anyARM']==0)).sum()),
          both_neg=int(((comp['amd_self']==0)&(comp['obj_anyARM']==0)).sum()))

# ---------- prespecified DCA + NNE (NHIS pooled, full) ----------
def dca_grid(y,lp,thr_list):
    y=np.asarray(y,float); p=1/(1+np.exp(-lp)); n=len(y); prev=y.mean(); rows=[]
    for pt in thr_list:
        flag=p>=pt; tp=((flag)&(y==1)).sum(); fp=((flag)&(y==0)).sum()
        nb=(tp-(fp*pt/(1-pt)))/n; allnb=prev-(1-prev)*pt/(1-pt)
        rows.append(dict(threshold=pt,nb_model=round(float(nb),5),nb_all=round(float(allnb),5),
            flagged_per_1000=int(round(1000*flag.mean(),0)),
            cases_identified_per_1000=int(round(1000*tp/n,1)),
            sensitivity=round(float(tp/y.sum()),3),
            cases_missed_per_1000=int(round(1000*(y.sum()-tp)/n,1))))
    return rows
dca_prespec=dca_grid(yext,lp_ext_all,[0.02,0.05,0.08,0.10])

J(dict(models=models,design=design,temporal=temporal,missingness=missingness,subgroups=subgroups,
       coef_table=coef_table,age_wald=age_wald,EPV=EPV,worked_example=worked,knots=kn,
       objective=dict(meta=obj,metrics=obj_metrics,concordance=conc),dca_prespec=dca_prespec),
   "harden_results.json")
print(json.dumps({"knots":kn,"nested":{k:models[k] for k in models if k not in('delta_auc_full_minus_ageonly',)},
                  "delta_auc":models['delta_auc_full_minus_ageonly'],
                  "temporal":temporal,"EPV":EPV,"age_wald":age_wald,"worked":worked,
                  "psu_ci":psu_ci,"obj":obj_metrics,"conc":conc},indent=2)[:4000])
