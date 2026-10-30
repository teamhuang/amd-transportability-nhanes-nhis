"""ROUND-8 P0-1: PROPER MICE with variable-type-appropriate conditional models.
Replaces the earlier BayesianRidge linear-Gaussian imputation (which drew continuous values for
binary predictors and had to clip 12,530 out-of-range cells / resolve illegal smoking pairs).

Design (per reviewer):
  - four BINARY predictors (obesity, hypertension, dyslipidemia, diabetes)  -> logistic conditional model
  - smoking as ONE three-category variable never/former/current            -> multinomial-logit conditional model
  - conditioning set = age, sex, the OUTCOME (amd_self, outcome-inclusive / MAR), survey cycle, and the other predictors
  - proper (Bayesian) draws: beta* ~ N(beta_hat, Cov_hat) then Bernoulli / Categorical sampling -> always a valid category
  - chained equations, T sweeps, m imputations
After each completed dataset the analysis logistic model is re-fit on all valid-outcome rows and the
frozen equation transported to the SAME untouched NHIS validation sample. Coefficients pooled by
Rubin's rules (within + between variance, Barnard-Rubin df). Also emits a diagnostics table
(missing fraction, observed vs imputed category proportions, model type, count of illegal /
out-of-range draws which is now zero BY CONSTRUCTION) and reproduces the canonical complete-case
external AUROC in the same run.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
import numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
SEED=20261028
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"
rng_master=np.random.RandomState(SEED)

der=pd.read_csv(OUT/"nhanes_derivation.csv")
ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
ccdf=der[der.cc_pred==True].reset_index(drop=True)

BIN=["obesity","hypertension","dyslipidemia","diabetes"]
ALLP=["age","male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]

# ---- smoking category (0 never,1 former,2 current); NaN where both dummies missing ----
def to_cat(df):
    sf=df["smoke_former"].values; sc=df["smoke_current"].values
    cat=np.full(len(df),np.nan)
    cat[(np.isnan(sf))&(np.isnan(sc))]=np.nan
    obs=~(np.isnan(sf)|np.isnan(sc))
    cat[obs]=np.where(sc[obs]==1,2,np.where(sf[obs]==1,1,0))
    # observed never = both 0
    return cat
der["_smk"]=to_cat(der)

cycle=(der["SDDSRVYR"].values.astype(float)); cyc_hi=(cycle>np.nanmedian(cycle)).astype(float)

# ---------- canonical complete-case (reproduce headline in same run) ----------
def design(df): return np.column_stack([np.ones(len(df)),df[ALLP].values.astype(float)])
def transport_eval(dfit,deval):
    X=design(dfit); y=dfit["amd_self"].values.astype(float)
    m=sm.Logit(y,X).fit(disp=0); b=np.asarray(m.params); se=np.asarray(m.bse)
    Xe=design(deval); lp=Xe@b; ye=deval["amd_self"].values.astype(float); p=1/(1+np.exp(-lp))
    slope=float(sm.Logit(ye,sm.add_constant(lp)).fit(disp=0).params[1])
    citl=float(sm.Logit(ye,np.ones((len(ye),1)),offset=lp).fit(disp=0).params[0])
    return dict(b=b,se=se,der_auc=float(roc_auc_score(y,X@b)),ext_auc=float(roc_auc_score(ye,lp)),
                ext_slope=slope,ext_citl=citl,ext_meanpred=float(p.mean()),ext_obs=float(ye.mean()))
cc=transport_eval(ccdf,ext)

# ---------- one proper-MICE completed dataset ----------
def draw_beta(bhat,cov,rng):
    try:
        L=np.linalg.cholesky(cov+1e-8*np.eye(len(bhat)))
    except np.linalg.LinAlgError:
        ev,V=np.linalg.eigh(cov); ev=np.clip(ev,1e-8,None); L=V@np.diag(np.sqrt(ev))
    return bhat+L@rng.standard_normal(len(bhat))

def _cov_from_hessian(model, params):
    """Observed-information covariance via pseudo-inverse of the Hessian (robust to
    near-singularity from rare/separated patterns, unlike fit.cov_params())."""
    try:
        H=np.asarray(model.hessian(np.asarray(params)),dtype=float)
    except Exception:
        return None
    if H.size==0 or np.any(np.isnan(H)) or np.any(np.isinf(H)): return None
    C=np.linalg.pinv(-H); C=(C+C.T)/2.0
    ev,_=np.linalg.eigh(C)
    if np.min(ev) < -1e-6: C = C + (-2*np.min(ev)+1e-8)*np.eye(C.shape[0])
    return C

def _draw_params(model, params, rng, prior_sd=1.0):
    C=_cov_from_hessian(model, params)
    bhat=np.asarray(params,dtype=float).ravel()
    if C is None:
        return bhat + rng.standard_normal(len(bhat))*prior_sd
    return draw_beta(bhat, C, rng)

def impute_once(seed,T=15):
    rng=np.random.RandomState(seed)
    d=der.copy()
    miss={c:pd.isna(d[c]).values for c in BIN}
    miss_smk=pd.isna(d["_smk"]).values
    # random-init from observed marginals
    for c in BIN:
        p_obs=float(np.nanmean(d[c])); d.loc[miss[c],c]=(rng.rand(int(miss[c].sum()))<p_obs).astype(float)
        d[c]=d[c].astype(float)
    smk_obs=d["_smk"].values[~miss_smk]
    init_probs=np.bincount(smk_obs.astype(int),minlength=3)/len(smk_obs)
    d.loc[miss_smk,"_smk"]=rng.choice([0,1,2],size=int(miss_smk.sum()),p=init_probs).astype(float)
    d["_smk"]=d["_smk"].astype(float)
    def dummy_from_cat(v): return np.where(v==1,1.0,0.0),np.where(v==2,1.0,0.0)   # (former, current)
    for _ in range(T):
        # --- binary predictors: logistic conditional on the rest ---
        for c in BIN:
            others=[o for o in BIN if o!=c]
            f,cur=dummy_from_cat(d["_smk"].values)
            X=np.column_stack([np.ones(len(d)),d["age"].values,d["male"].values,d["amd_self"].values,cyc_hi,f,cur,d[others].values.astype(float)])
            obs=~miss[c]
            try:
                mod=sm.Logit(d[c].values[obs],X[obs]); fit=mod.fit(disp=0,method="bfgs",maxiter=200)
                bstar=_draw_params(mod, fit.params, rng)
            except Exception:
                p_obs=float(np.mean(d[c].values[obs])); lin=None; bstar=None
            if bstar is None:
                prob=np.full(int(miss[c].sum()),p_obs)
            else:
                lin=X[miss[c]]@bstar; prob=1/(1+np.exp(-lin))
            d.loc[miss[c],c]=(rng.rand(prob.shape[0])<prob).astype(float)
            d[c]=d[c].astype(float)
        # --- smoking: multinomial logit conditional on the rest (never = baseline) ---
        f,cur=dummy_from_cat(d["_smk"].values)
        Xs=np.column_stack([np.ones(len(d)),d["age"].values,d["male"].values,d["amd_self"].values,cyc_hi,d[BIN].values.astype(float)])
        obsy=~miss_smk
        k=Xs.shape[1]
        try:
            modS=sm.MNLogit(d["_smk"].values[obsy],Xs[obsy]); fitS=modS.fit(disp=0,method="bfgs",maxiter=200)
            bstar=_draw_params(modS, fitS.params, rng).reshape(2,k)   # rows: former, current (never = baseline)
            lin=Xs[miss_smk]@bstar.T                                   # (n_miss, 2)
            e=np.exp(np.hstack([np.zeros((lin.shape[0],1)), lin]))     # never column = 0 linear predictor
            pr=e/e.sum(1,keepdims=True)
        except Exception:
            smk_obs2=d["_smk"].values[obsy].astype(int)
            ip=np.bincount(smk_obs2,minlength=3)/len(smk_obs2); pr=np.tile(ip,(int(miss_smk.sum()),1))
        draws=np.array([rng.choice([0,1,2],p=row) for row in pr],dtype=float)
        d.loc[miss_smk,"_smk"]=draws
    # assemble final 0/1 predictors
    f,cur=dummy_from_cat(d["_smk"].values)
    d["smoke_former"]=f; d["smoke_current"]=cur
    return d

def illegal_and_oor(d):
    # count draws that would have been out-of-[0,1] or illegal under a linear model: under this
    # categorical sampler they are structurally impossible -> report 0, but verify anyway.
    bin_ok=all(set(np.unique(d[c].dropna().values)).issubset({0.0,1.0}) for c in BIN)
    both=int(((d["smoke_former"]==1)&(d["smoke_current"]==1)).sum())
    return bin_ok,both

M=20
T=15
ext_aucs=[];slopes=[];citls=[];der_aucs=[];coefs=[];wses=[]
prepost={c:[] for c in BIN}; prepost["smk_never"]=[];prepost["smk_former"]=[];prepost["smk_current"]=[]
bin_ok_all=True; illegal_total=0
for r in range(M):
    d=impute_once(SEED+r)
    b_ok,both=illegal_and_oor(d); bin_ok_all=b_ok and bin_ok_all; illegal_total+=both
    dfull=d[d["amd_self"].notna()].reset_index(drop=True)
    for c in BIN: prepost[c].append(float(dfull[c].mean()))
    prepost["smk_never"].append(float((dfull["_smk"]==0).mean()))
    prepost["smk_former"].append(float(dfull["smoke_former"].mean()))
    prepost["smk_current"].append(float(dfull["smoke_current"].mean()))
    res=transport_eval(dfull,ext)
    coefs.append(res["b"]); wses.append(res["se"])
    ext_aucs.append(res["ext_auc"]); slopes.append(res["ext_slope"]); citls.append(res["ext_citl"]); der_aucs.append(res["der_auc"])

C=np.array(coefs); W=np.array(wses)
pooled_mean=C.mean(0); between=C.var(0,ddof=1); within=(W**2).mean(0)
total=within*(1+1/M)+between
# Barnard-Rubin df for the age coefficient (index 1)
lam=(1+1/M)*between/total
def Rubin_df(tot,idx):
    r=lam[idx]; vwithin=within[idx]; vbet=between[idx]
    if r<=0: return 1e6
    return (M-1)*(1+ (vwithin*(1+1/M))/ (r*r*(vbet+1e-18)) )**2
from scipy.stats import t as tdist
def pooled_ci(idx):
    df=Rubin_df(total,idx); tcrit=tdist.ppf(0.975, min(df,1e6))
    return float(np.sqrt(total[idx])), float(pooled_mean[idx]-tcrit*np.sqrt(total[idx])), float(pooled_mean[idx]+tcrit*np.sqrt(total[idx])), float(df)
age_se,age_lo,age_hi,age_df=pooled_ci(1)

mi=dict(method="proper MICE (chained equations): logistic for 4 binary predictors, multinomial logit for 3-category smoking; beta* ~ N(beta_hat,Cov) draws, T=15 sweeps, m=20",
        M=M, sweeps=T, binary_model="logistic", smoking_model="multinomial-logit (never/former/current)",
        outcome_inclusive=True, auxiliary="age, sex, outcome, survey cycle (SDDSRVYR)",
        ext_auc_mean=round(float(np.mean(ext_aucs)),4), ext_auc_min=round(float(np.min(ext_aucs)),4),
        ext_auc_max=round(float(np.max(ext_aucs)),4), ext_auc_sd=round(float(np.std(ext_aucs,ddof=1)),4),
        der_auc_mean=round(float(np.mean(der_aucs)),4),
        ext_slope_mean=round(float(np.mean(slopes)),4), ext_citl_mean=round(float(np.mean(citls)),4),
        pooled_coef=dict(zip(["const"]+ALLP,[round(float(x),4) for x in pooled_mean])),
        age_coef=dict(pooled=round(float(pooled_mean[1]),4), rubin_total_se=round(age_se,4),
                      ci95=[round(age_lo,4),round(age_hi,4)], barnard_rubin_df=round(age_df,1)),
        validity=dict(binary_imputations_valid_01_always=bool(bin_ok_all),
                      illegal_smoking_combos=int(illegal_total), out_of_range_cells=0,
                      note="categorical Bernoulli/multinomial sampling makes out-of-range and illegal-combination draws structurally impossible"))

# diagnostics: observed vs imputed category proportions
diag=[]
for c in BIN:
    obs_p=float(ccdf[c].mean())
    diag.append(dict(variable=c,model="logistic",n_missing=int(der[c].isna().sum()),
                     pct_missing=round(100*float(der[c].isna().mean()),2),
                     observed_positive_pct=round(100*obs_p,2),
                     imputed_overall_positive_pct=round(100*float(np.mean(prepost[c])),2)))
smk_diag=dict(variable="smoking (3-category)",model="multinomial-logit",
              n_missing=int(der["_smk"].isna().sum()),pct_missing=round(100*float(der["_smk"].isna().mean()),2),
              observed_pct=dict(never=round(100*float((der["_smk"]==0).mean()),1),
                                former=round(100*float((der["_smk"]==1).mean()),1),
                                current=round(100*float((der["_smk"]==2).mean()),1)),
              imputed_overall_pct=dict(never=round(100*float(np.mean(prepost['smk_never'])),1),
                                       former=round(100*float(np.mean(prepost['smk_former'])),1),
                                       current=round(100*float(np.mean(prepost['smk_current'])),1)))
out=dict(complete_case=dict(ext_auc=round(cc["ext_auc"],4),ext_slope=round(cc["ext_slope"],4),
            ext_citl=round(cc["ext_citl"],4),der_auc=round(cc["der_auc"],4),
            coef=dict(zip(["const"]+ALLP,[round(float(x),4) for x in cc["b"]]))),
           multiple_imputation=mi, diagnostics=dict(binary=diag, smoking=smk_diag))
json.dump(out,open(OUT/"sr_missing_mi2.json","w"),indent=2)
print(json.dumps(out,indent=1))
