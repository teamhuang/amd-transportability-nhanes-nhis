"""Topic 1 main analysis: derivation (NHANES) + internal validation + external validation
(NHIS pooled / per-year / per-region) + objective-grading sensitivity. Frozen admission_record.md.

Unweighted logistic derivation (TRIPOD individual-risk convention); design weights carried for
weighted descriptives/sensitivity. Seeds fixed; BLAS threads locked before numpy.
Outputs: model_coefficients.csv, internal_validation.json, external_metrics.csv, sensitivity_objective.json, findings.json
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[v]="1"
import numpy as np, pandas as pd, json, tempfile
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit
from statsmodels.stats.multitest import multipletests
import patsy
from sklearn.metrics import roc_auc_score, brier_score_loss

SEED=20261028
OUT=r"H:\眼科论文\题目4\analysis\output"
PRED=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
RHS="cr(age, df=3) + "+" + ".join(PRED)
YCOL="amd_self"

def atomic_json(obj,path):
    tmp=tempfile.NamedTemporaryFile('w',delete=False,suffix='.json',dir=OUT,newline='').name
    with open(tmp,'w') as f: json.dump(obj,f,indent=2)
    os.replace(tmp,path)

def wmean(s,w):
    w=w/np.nansum(w); return float(np.nansum(s*w))

def auc_boot(y,p,seed=SEED,nboot=1000):
    auc=roc_auc_score(y,p); rng=np.random.RandomState(seed); n=len(y); boots=[]
    for _ in range(nboot):
        idx=rng.randint(0,n,n)
        yy=y[idx]; 
        if yy.sum()==0 or yy.sum()==len(yy): continue
        try: boots.append(roc_auc_score(yy,p[idx]))
        except: pass
    lo,hi=np.percentile(boots,[2.5,97.5])
    return dict(auc=round(float(auc),4),ci=[round(float(lo),4),round(float(hi),4)],nboot=len(boots))

def calib(y,lp):
    # slope: logit(y)~lp ; intercept-in-the-large: logit(y)~1 offset=lp
    import numpy as _np
    X=sm.add_constant(lp)
    fit=Logit(_np.asarray(y,dtype=float),X).fit(disp=0)
    slope=float(fit.params[1]); slope_ci=[float(x) for x in fit.conf_int()[1]]
    Xint=_np.ones((len(y),1))
    fit2=Logit(_np.asarray(y,dtype=float),Xint,offset=lp).fit(disp=0)
    c_intercept=float(fit2.params[0]); ci=[float(x) for x in fit2.conf_int()[0]]
    return dict(slope=round(slope,4),slope_ci=[round(s,4) for s in slope_ci],
               calib_intercept=round(c_intercept,4),intercept_ci=[round(float(x),4) for x in ci])

def net_benefit(y,p,pt):
    y=np.asarray(y,dtype=float); p=np.asarray(p,dtype=float); n=len(y); prev=y.mean()
    tp=np.sum((p>=pt)&(y==1)); fp=np.sum((p>=pt)&(y==0))
    nb=(tp-(fp*pt/(1-pt)))/n
    all_nb=prev-(1-prev)*pt/(1-pt)
    return float(nb),float(all_nb)

def dca(y,lp_probs,thresholds):
    """net benefit of model vs all vs none at each threshold."""
    y=np.asarray(y,dtype=float); p=np.asarray(lp_probs,dtype=float); n=len(y); prev=y.mean()
    rows=[]
    for pt in thresholds:
        if pt<=0 or pt>=1: continue
        tp=np.sum((p>=pt)&(y==1)); fp=np.sum((p>=pt)&(y==0))
        nb=(tp - (fp*pt/(1-pt)))/n
        all_nb=prev-(1-prev)*pt/(1-pt)
        rows.append(dict(threshold=round(pt,3),net_benefit_model=round(float(nb),5),
                         net_benefit_all=round(float(all_nb),5),net_benefit_none=0.0))
    return rows

# ---------------- LOAD ----------------
der=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv"))
der=der[der['cc_pred']==True].copy()
ext=pd.read_csv(os.path.join(OUT,"nhis_external.csv"))
ext=ext[ext['cc_pred']==True].copy()

# ---------------- DERIVATION ----------------
yder,Xder=patsy.dmatrices(f"{YCOL} ~ {RHS}", der, return_type="dataframe")
design_info=Xder.design_info
model=Logit(yder.values.ravel(), Xder.values).fit(disp=0)
params=Xder.columns.tolist(); coef=model.params
se=model.bse; pvals=model.pvalues; conf=model.conf_int()
n=len(yder); yhat=model.predict(Xder.values)
print("DERIVATION converged:",bool(model.mle_retvals['converged']),"iter",int(model.mle_retvals.get('iterations',-1)))

# coefficients table (OR) on predictor terms
rows=[]
for i,name in enumerate(params):
    rows.append(dict(term=name,beta=float(coef[i]),se=float(se[i]),
        OR=float(np.exp(coef[i])),OR_lo=float(np.exp(conf[i][0])),OR_hi=float(np.exp(conf[i][1])),
        p_raw=float(pvals[i])))
coefdf=pd.DataFrame(rows)
# Holm across predictor terms (exclude Intercept)
maskp=coefdf['term']!='Intercept'
rej,holm,_,_=multipletests(coefdf.loc[maskp,'p_raw'].values,method='holm')
coefdf.loc[maskp,'p_holm']=holm
coefdf.to_csv(os.path.join(OUT,"model_coefficients.csv"),index=False)

# global age effect test: joint Wald on the cr() terms? report each spline term p.
auc_apparent=roc_auc_score(yder.values.ravel(),yhat)

# ---------------- INTERNAL VALIDATION (bootstrap optimism, Harrell) ----------------
rng=np.random.RandomState(SEED); B=500
opt_auc=[]; opt_slope=[]
Xv=Xder.values; Yv=yder.values.ravel().astype(float)
for b in range(B):
    idx=rng.randint(0,n,n)
    try:
        mb=Logit(Yv[idx],Xv[idx]).fit(disp=0)
        lp_b=mb.predict(Xv)              # test on original sample
        lp_bb=mb.predict(Xv[idx])        # apparent on bootstrap
        a_test=roc_auc_score(Yv,lp_b); a_app=roc_auc_score(Yv[idx],lp_bb)
        opt_auc.append(a_app-a_test)
        # calibration slope on original using bootstrap model
        sl=Logit(Yv.astype(float),sm.add_constant(lp_b)).fit(disp=0).params[1]
        opt_slope.append(1.0-sl)
    except Exception:
        continue
auc_corrected=auc_apparent-np.mean(opt_auc)
slope_apparent=1.0
calib_app=calib(Yv.astype(float), np.log(np.clip(yhat,1e-9,1-1e-9)/(1-np.clip(yhat,1e-9,1-1e-9))))
# bootstrap calibration-slope optimism is unstable (degenerate resamples); use MEDIAN optimism
opt_slope_med=float(np.median(opt_slope)) if opt_slope else 0.0
slope_corrected=slope_apparent-opt_slope_med
internal=dict(apparent_auc=round(float(auc_apparent),4),
              optimism_auc=round(float(np.mean(opt_auc)),4),
              corrected_auc=round(float(auc_corrected),4),
              n_boot=len(opt_auc), B=B, seed=SEED,
              note="internal calibration slope optimism is unstable at 214 events/10 terms and is intentionally NOT reported; calibration is assessed in the external cohorts (Van Houwelingen slope + intercept-in-the-large).",
              derivation_calib=calib_app,
              brier=round(float(brier_score_loss(Yv,yhat)),5))
atomic_json(internal,os.path.join(OUT,"internal_validation.json"))

# ---------------- EXTERNAL VALIDATION ----------------
# frozen LP: apply derivation design_info to each validation subset
def lp_for(subset):
    Xnew=patsy.build_design_matrices([design_info], subset, return_type="dataframe")[0]
    return Xnew.values@coef, np.array(Xnew.columns.tolist())

ext_sets={}
ext_sets["NHIS_pooled"]=ext
for y in (2016,2017): ext_sets[f"NHIS_{y}"]=ext[ext['survey_year']==y]
region_names={1:"Northeast",2:"Midwest",3:"South",4:"West"}
for r,nm in region_names.items():
    sub=ext[(ext['region']==r)]
    if len(sub)>500: ext_sets[f"NHIS_region_{nm}"]=sub

metrics=[]
for name,sub in ext_sets.items():
    sub=sub.copy()
    lp,_=lp_for(sub)
    y=sub[YCOL].values.astype(float)
    p=1/(1+np.exp(-lp))
    a=auc_boot(y,lp)
    c=calib(y,lp)
    prev=float(y.mean())
    thresholds=sorted(set([0.01,0.02,0.03,prev,round(prev*0.5,4),round(prev*1.5,4)]))
    d=dca(y,p,thresholds)
    nb_prev,all_prev=net_benefit(y,p,prev)
    metrics.append(dict(cohort=name,n=int(len(sub)),events=int(y.sum()),
        prevalence=round(prev,5),
        nb_at_prevalence=round(nb_prev,5), nb_all_at_prevalence=round(all_prev,5),
        auc=a['auc'],auc_lo=a['ci'][0],auc_hi=a['ci'][1],
        calib_slope=c['slope'],calib_slope_lo=c['slope_ci'][0],calib_slope_hi=c['slope_ci'][1],
        calib_intercept=c['calib_intercept'],intercept_lo=c['intercept_ci'][0],intercept_hi=c['intercept_ci'][1],
        brier=round(float(brier_score_loss(y,p)),5),
        mean_pred_prob=round(float(p.mean()),5),
        dca=d))
mrows=[]
for m in metrics:
    m2={k:v for k,v in m.items() if k!='dca'}
    # net benefit at prevalence threshold
    for r in m['dca']:
        if abs(r['threshold']-m['prevalence'])<1e-9:
            m2['nb_at_prevalence']=r['net_benefit_model']; m2['nb_all_at_prevalence']=r['net_benefit_all']
    mrows.append(m2)
pd.DataFrame(mrows).to_csv(os.path.join(OUT,"external_metrics_flat.csv"),index=False)
with open(os.path.join(OUT,"external_metrics_full.json"),"w") as f: json.dump(metrics,f,indent=2)

# full-recalibration comparison for NHIS_pooled (intercept-only and intercept+slope)
sub=ext_sets["NHIS_pooled"]; lp,_=lp_for(sub); y=sub[YCOL].values.astype(float)
auc_pooled=auc_boot(y,lp)

# ---- dump plotting-ready linear predictors (frozen model) for each external set ----
dump=[]
for name,sb in ext_sets.items():
    l,_x=lp_for(sb); pp=1/(1+np.exp(-l))
    t=pd.DataFrame(dict(cohort=name,lp=l,pred=pp,y=sb[YCOL].values.astype(float)))
    if 'region' in sb: t['region']=sb['region'].values
    if 'survey_year' in sb: t['survey_year']=sb['survey_year'].values
    dump.append(t)
pd.concat(dump,ignore_index=True).to_csv(os.path.join(OUT,"external_lp.csv"),index=False)
# derivation predicted probabilities + observed for calibration plot
dcal=pd.DataFrame(dict(lp=Xder.values@coef,pred=yhat,y=Yv,weight=der['WTMEC2YR'].values))
dcal.to_csv(os.path.join(OUT,"derivation_pred.csv"),index=False)

# ---- age dose-response: OR vs age (ref=55) holding others at sample means ----
from patsy import dmatrix
mean_cols={p:float(np.mean(der[p])) for p in PRED}
ages=list(range(45,86,5)); age_or=[]
# reconstruct design column order to build a new row
def row_for(age):
    basis=dmatrix("cr(x, df=3)", {"x":np.array([age])}, design_info=None)  # not used; handled below
    return None
# simpler: evaluate predicted prob on a synthetic person with mean covariates at each age using the SAME design_info
base=pd.DataFrame({c:[mean_cols[c] for _ in ages] for c in PRED})
base['age']=ages
for c in ['amd_self']: base[c]=0
Xa=patsy.build_design_matrices([design_info], base, return_type="dataframe")[0]
pa=1/(1+np.exp(-(Xa.values@coef)))
ref_idx=list(ages).index(55) if 55 in ages else 2
odds=pa/(1-pa)
age_or=dict(ages=ages, prob=[round(float(v),5) for v in pa], OR_vs_ref55=[round(float(odds[i]/odds[ref_idx]),3) for i in range(len(ages))])
atomic_json(age_or,os.path.join(OUT,"age_dose.json"))
# simple integer points score (per-binary-predictor, ref = 55y female never-smoker normal everything)
base_ref=pd.DataFrame({c:[0.0]*1 for c in PRED}); base_ref['age']=55.0
Xr=patsy.build_design_matrices([design_info], base_ref, return_type="dataframe")[0].values[0]
lp0=float(Xr@coef)
pts=dict(intercept_lp_at_ref=round(lp0,4),
         points={p:int(round(coef[list(params).index(p)])) for p in PRED},
         betas={p:round(float(coef[list(params).index(p)]),4) for p in PRED})
atomic_json(pts,os.path.join(OUT,"score_points.json"))

# ---------------- OBJECTIVE-GRADING SENSITIVITY (NHANES) ----------------
sens=dict()
for ey,lab in [('obj_anyARM','objective_anyARM'),('obj_lateARM','objective_lateARM')]:
    ssub=der[der[ey].notna()].copy()
    Xs=patsy.build_design_matrices([design_info], ssub, return_type="dataframe")[0]
    lp_s=Xs.values@coef
    yv=ssub[ey].values.astype(float)
    if np.nansum(yv)>10 and not np.isnan(yv).any():
        a=auc_boot(yv,lp_s)
        sens[lab]=dict(n=int(len(yv)),events=int(yv.sum()),prev=round(float(yv.mean()),5),
                       auc_frozen_selfreport_score_on=a['auc'],ci=a['ci'])
# concordance self-report vs objective anyARM
comp=der[[YCOL,'obj_anyARM']].dropna()
sens['concordance']={'n':int(len(comp)),'both_yes':int(((comp[YCOL]==1)&(comp['obj_anyARM']==1)).sum()),
    'self_only':int(((comp[YCOL]==1)&(comp['obj_anyARM']==0)).sum()),
    'obj_only':int(((comp[YCOL]==0)&(comp['obj_anyARM']==1)).sum()),
    'neither':int(((comp[YCOL]==0)&(comp['obj_anyARM']==0)).sum())}
# re-fit score on objective endpoint (anyARM)
try:
    der_obj=der[der['obj_anyARM'].notna()].copy()
    yo2,Xo2=patsy.dmatrices("obj_anyARM ~ "+" + ".join(["cr(age,df=3)"]+PRED), der_obj, return_type="dataframe")
    mo=Logit(yo2.values.ravel(),Xo2.values).fit(disp=0)
    yo=yo2.values.ravel(); pmo=mo.predict(Xo2.values)
    sens['obj_model_auc']=round(float(roc_auc_score(yo,pmo)),4)
    sens['obj_model_n']=int(len(yo)); sens['obj_model_events']=int(yo.sum())
except Exception as e:
    sens['obj_model_auc']=f"err:{e}"
atomic_json(sens,os.path.join(OUT,"sensitivity_objective.json"))

# ---------------- FINDINGS SUMMARY ----------------
findings=dict(
  derivation_n=int(n), events=int(yder.values.sum()), prev=round(float(yder.values.mean()),5),
  EPV_per_predictor=round(float(yder.values.sum()/len(PRED)),1),
  apparent_auc=internal['apparent_auc'], corrected_auc=internal['corrected_auc'],
  external=[dict(cohort=m['cohort'],n=m['n'],prev=m['prevalence'],auc=m['auc'],
                 auc_lo=m['auc_lo'],auc_hi=m['auc_hi'],
                 calib_slope=m['calib_slope'],calib_intercept=m['calib_intercept'],
                 nb_at_prev=m['nb_at_prevalence']) for m in metrics],
  objective_sensitivity=sens)
atomic_json(findings,os.path.join(OUT,"findings.json"))
print("=== DERIVATION coefficients (OR) ===")
print(coefdf[['term','OR','OR_lo','OR_hi','p_raw','p_holm']].to_string(index=False))
print("\n=== INTERNAL ===")
print(json.dumps(internal,indent=2))
print("\n=== EXTERNAL ===")
print(pd.DataFrame(mrows)[['cohort','n','events','prevalence','auc','auc_lo','auc_hi','calib_slope','calib_intercept','brier','nb_at_prevalence']].to_string(index=False))
print("\n=== OBJECTIVE SENSITIVITY ===")
print(json.dumps(sens,indent=2))
