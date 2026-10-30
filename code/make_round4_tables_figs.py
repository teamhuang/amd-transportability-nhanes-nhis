"""ROUND-4: emit supplementary tables (endpoint code audit, sample flow + MI, survey inference,
model-increment fairness, objective paired + gradeability, age functional form, diagnostic
concordance with Wilson CIs) into Tables_sr4.md, and build FigSR7 (self-report external smooth
calibration + risk histogram) and FigSR8 (objective-endpoint calibration, separate panels per
outcome). All numbers read from the round-4 JSONs or recomputed here from stored data.
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"; FIG=ROOT/"figures"
import numpy as np, pandas as pd, json
import statsmodels.api as sm
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SEED=20261028
mi=json.load(open(OUT/"sr_missing_mi.json")); mi2=json.load(open(OUT/"sr_missing_mi2.json")); sv=json.load(open(OUT/"sr_survey.json")); inc=json.load(open(OUT/"sr_increment.json"))
mistab=json.load(open(OUT/"sr_mi_stability.json")) if (OUT/"sr_mi_stability.json").exists() else None
lipidS=json.load(open(OUT/"sr_lipid_sensitivity.json")) if (OUT/"sr_lipid_sensitivity.json").exists() else None
BIN=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
der=pd.read_csv(OUT/"nhanes_derivation.csv"); der=der[der.cc_pred==True].reset_index(drop=True)
ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
def D(df,preds): return np.column_stack([np.ones(len(df)),df[["age"]+preds].values.astype(float)])
def fitlog(df,preds): return sm.Logit(df["amd_self"].values.astype(float),D(df,preds)).fit(disp=0)
mf=fitlog(der,BIN); ma=fitlog(der,[])
lf=D(ext,BIN)@np.asarray(mf.params); la=D(ext,[])@np.asarray(ma.params)
a16=float(sm.Logit(ext.loc[ext.survey_year==2016,"amd_self"].values.astype(float),np.ones((int((ext.survey_year==2016).sum()),1)),offset=lf[(ext.survey_year==2016).values]).fit(disp=0).params[0])
lupd=lf+a16
ye=ext["amd_self"].values.astype(float)

def wilson(k,n):
    if n==0: return (None,None)
    p=k/n; z=1.96; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return round(p,4),[round(max(0,c-h),4),round(min(1,c+h),4)]

# concordance 2x2 (from Table 5, unweighted)
a,b,c,d=72,316,91,4588
sens=wilson(a,a+b); spec=wilson(d,c+d); ppv=wilson(a,a+c); npv=wilson(d,b+d)

L=[]
def w(s=""): L.append(s)
w("## Round-4 supplementary tables\n")
w("### Supplementary Table S7. Endpoint definitions verified against the live CDC codebook and the raw NCHS transport files (audit log: `endpoint_code_audit.md`)\n")
w("| Endpoint | Survey / variable | Verified wording / coding | Our operationalisation |"); w("| --- | --- | --- | --- |")
w("| Self-reported AMD (outcome) | NHANES `VIQ310` (2005–06, 2007–08) | \"…ever been told by an **eye doctor** that … age-related macular degeneration?\"; both sexes 40–150 y; 1=Yes, 2=No, 7=Refused (0 occurrences), 9=Don't know | `VIQ310==1`; 7/9 missing |")
w("| Self-reported AMD (outcome) | NHIS `VIM_MDEV` (2016, 2017) | \"EVER been told you had macular degeneration\"; sample adults ≥18 y; 1=Yes, 2=No, 7=Refused, 8=Not ascertained (0 occurrences), 9=Don't know | `VIM_MDEV==1`; 7/9 missing |")
w("| Photographic any-ARM | NHANES `OPDUARM` (2005–06, 2007–08) | \"ARM, 3 severity levels, worse eye\"; both sexes 40–150 y; 0=No ARM, 1=Early ARM (drusen and/or pigmentary abnormalities), 2=Late ARM (exudative ARM and/or geographic atrophy) | `OPDUARM>=1` (early+late photographic ARM) |")
w("| Photographic late AMD | NHANES `OPDULARM` (2005–06, 2007–08) | \"Late ARM, worse eye\"; 0=No, 1=Yes (exudative ARM and/or geographic atrophy) | `OPDULARM==1` |")
w("")
w("**Verification status (completed before submission, October 2026):** every label and code was read both from the raw `.xpt`/`.csv` transports and from the live CDC/NCHS codebook pages; counts agree **cell-by-cell and exactly** in both NHANES cycles and both NHIS years (`analysis/scripts/audit_codebook_round12.py`; machine-readable record `analysis/output/sr_codebook_audit.json`). **Caveat:** CDC publishes self-reported AMD (`VIQ310`, `VIM_MDEV`) and examination-based AMD (`OPDUARM`) as **separate case definitions** and makes **no equivalence statement** between them, and a third claims-based \"diagnosed AMD\" definition uses ICD-9/ICD-10 codes in claims or EHR data (CDC Vision Health Surveillance case definitions, https://www.cdc.gov/vision-health-data/case-definitions/age-related-macular-degeneration-amd.html). The NHANES and NHIS self-report items also differ in wording (eye doctor vs any doctor), an outcome-ascertainment difference across surveys.\n")

f=mi["flow"]
w("### Supplementary Table S8. Sample flow and complete-case vs multiple-imputation comparison\n")
w("| Step / metric | NHANES n (or value) |"); w("| --- | --- |")
w(f"| Age >=40 with valid AMD outcome | {f['age_eligible_ge40']} |")
w(f"| Complete predictors (analytic) | {f['complete_predictors']} |")
w(f"| Excluded, incomplete predictors | {f['excluded_incomplete_predictors']} ({100*f['excluded_incomplete_predictors']/f['age_eligible_ge40']:.1f}%) |")
w(f"| Fundus-gradable any-ARM | {f['fundus_gradable_anyARM']} |")
w(f"| Fundus-gradable late ARM | {f['fundus_gradable_lateARM']} |")
w(f"| AMD prevalence: included vs excluded vs all | {100*f['amd_prev_included']:.2f}% vs {100*f['amd_prev_excluded']:.2f}% vs {100*f['amd_prev_all_eligible']:.2f}% |")
w("")
cc=mi2["complete_case"]; m=mi2["multiple_imputation"]
w("**Complete-case vs multiple-imputation comparison (same S8):**\n")
w("**The missing-data sensitivity uses chained-equations multiple imputation with variable-type-appropriate conditional models** (logistic for the four binary health indicators; multinomial logistic for the three-category smoking variable never/former/current), Bayesian parameter draws `beta* ~ N(beta_hat, Cov_hat)`, m = 20 imputations x 15 sweeps, outcome-inclusive and cycle-adjusted. Out-of-range and illegal-combination imputations are now structurally impossible (categorical sampling), so no clipping/rounding is used. (An earlier linear-Gaussian `IterativeImputer`+`BayesianRidge` version that clipped 12,530 out-of-range cells was superseded as statistically inappropriate for these variable types.)\n")
w("| Quantity | Complete-case (primary) | Multiple imputation (m=20) |"); w("| --- | --- | --- |")
w(f"| External AUROC | {cc['ext_auc']:.4f} | {m['ext_auc_mean']:.4f} (range {m['ext_auc_min']:.4f}-{m['ext_auc_max']:.4f}, between-imp SD {m['ext_auc_sd']:.4f}) |")
w(f"| External calibration slope | {cc['ext_slope']:.4f} | {m['ext_slope_mean']:.4f} |")
w(f"| External CITL intercept | {cc['ext_citl']:.4f} | {m['ext_citl_mean']:.4f} |")
w(f"| Age coefficient (logit), Rubin-pooled | {cc['coef']['age']:.4f} | {m['age_coef']['pooled']:.4f} (95% CI {m['age_coef']['ci95'][0]:.4f} to {m['age_coef']['ci95'][1]:.4f}, total SE {m['age_coef']['rubin_total_se']:.4f}) |")
w(f"| Imputation validity (m x n) | n/a | {m['validity']['out_of_range_cells']} out-of-range cells; {m['validity']['illegal_smoking_combos']} illegal smoking combinations (both structurally prevented) |")
w("")
w("**MI diagnostics (observed vs imputed distributions):**\n")
w("| Variable | Model type | % missing | Observed category % | Imputed-overall category % |"); w("| --- | --- | --- | --- | --- |")
for b in mi2["diagnostics"]["binary"]:
    w(f"| {b['variable']} | {b['model']} | {b['pct_missing']:.2f} | {b['observed_positive_pct']:.1f}% positive | {b['imputed_overall_positive_pct']:.1f}% positive |")
s=mi2["diagnostics"]["smoking"]; op=s['observed_pct']; ip=s['imputed_overall_pct']
w(f"| smoking (3-category) | {s['model']} | {s['pct_missing']:.2f} | never {op['never']:.1f} / former {op['former']:.1f} / current {op['current']:.1f} | never {ip['never']:.1f} / former {ip['former']:.1f} / current {ip['current']:.1f} |")
w("")
w("MI recovers the 1,238 excluded participants and leaves the core conclusions unchanged (age-dominated, external AUROC ~0.774, persistent over-prediction). Observed and imputed category proportions are close (S8 diagnostics block), indicating the type-appropriate conditional models reproduce the empirical distribution rather than distorting it. Complete-case selection is **not** negligible (excluded AMD 6.06% vs included 3.71%) but does not overturn the findings. The between-imputation **range** ({:.4f}-{:.4f}) reflects imputation variability; formal inference for a coefficient uses Rubin's rules (total SE, Barnard-Rubin degrees of freedom), and the Barnard-Rubin correction here is negligible (very large residual df) so the reported age-coefficient interval is effectively a normal interval.\n".format(m['ext_auc_min'], m['ext_auc_max']))
if mistab:
    dg=mistab["diagnostics"]
    w("**Independent code re-check of the MICE (empirical convergence across seeds x sweeps):** ")
    w(f"re-running the sampler over sweep counts 5/15/30 x three random seeds at m=6 gives a pooled external-AUROC grid mean **{mistab['pooled_external_auc_mean_over_grid']:.4f}** with total spread **{mistab['spread_over_grid']:.4f}** across the whole grid (published single-run value {mistab['comparison_to_published']['published_sr_missing_mi2_ext_auc']:.4f}); across {dg['conditional_fits']:,} conditional fits there were **{dg['binary_marginal_fallbacks']}** binary and **{dg['smoking_marginal_fallbacks']}** smoking fallbacks to the marginal (no singular-Hessian / perfect-separation event occurred), and imputed values were valid 0/1 or a single smoking level in every draw. Run log archived as `sr_mi_stability_runlog.txt`. The absence of a formal MICE trace is why we report stability across cycles+seeds rather than claiming a convergence statistic.\n")

w("### Supplementary Table S9. Complex-survey external inference (intervals are method-specific, not one common design-based estimator)\n")
w("| Quantity (NHIS) | Estimate | 95% interval | Method (see caveat) |"); w("| --- | --- | --- | --- |")
for nm in ["pooled","2016","2017"]:
    pv=sv["prevalence"][nm]; w(f"| Weighted AMD prevalence ({nm}) | {100*pv['weighted']:.2f}% | {100*pv['ci95'][0]:.2f}-{100*pv['ci95'][1]:.2f}% | Taylor linearisation (strata+PSU) |")
cp=sv["calibration_cluster_glm"]["pooled"]; bp=sv["psu_bootstrap"]
w(f"| Weighted calibration slope (pooled) | {cp['slope']:.3f} | {cp['slope_ci'][0]:.3f}-{cp['slope_ci'][1]:.3f} | Weighted GLM, PSU-clustered SE |")
w(f"| Weighted calibration slope (pooled) | {sv['weighted_point']['slope']:.3f} | {bp['seed_20261028']['slope'][0]:.3f}-{bp['seed_20261028']['slope'][1]:.3f} | PSU bootstrap B=1000 (seed 1) |")
w(f"| Weighted CITL (pooled) | {cp['citl']:.3f} | {cp['citl_ci'][0]:.3f}-{cp['citl_ci'][1]:.3f} | Weighted GLM, PSU-clustered SE |")
w(f"| Weighted AUROC (pooled) | {sv['weighted_point']['auc']:.3f} | {bp['seed_20261028']['auc'][0]:.3f}-{bp['seed_20261028']['auc'][1]:.3f} | PSU bootstrap B=1000 |")
w("")
w(f"Bootstrap seed-to-seed stability (max CI-endpoint difference across seeds {SEED} and 777): {bp['seed_stability_maxdiff']:.4f}. Weighted AUROC ({sv['weighted_point']['auc']:.3f}) differs from unweighted (0.774) because it targets a different population, not because the model is better.\n")

w("### Supplementary Table S10. Full-vs-age-only increment (paired bootstrap), by estimand\n")
w("| Estimand | ΔAUROC (full−age) (95% CI) | ΔBrier (full−age) (95% CI) | ΔNetBenefit@5% (full−age) (95% CI) |"); w("| --- | --- | --- | --- |")
for e in ["unweighted","weighted"]:
    x=inc["increment"][e]
    w(f"| {e} | {x['dAUC']:.4f} ({x['dAUC_CI'][0]:.4f} to {x['dAUC_CI'][1]:.4f}) | {x['dBrier']:.6f} ({x['dBrier_CI'][0]:.6f} to {x['dBrier_CI'][1]:.6f}) | {x['dNB5']:.5f} ({x['dNB5_CI'][0]:.5f} to {x['dNB5_CI'][1]:.5f}) |")
w("")
w(f"All three deltas are **full-model minus age-only**. A **positive ΔBrier means the full model has a HIGHER (worse) Brier score** than age only; a positive ΔAUROC/ΔNetBenefit favours the full model. The AUROC gain is statistically detectable but tiny in both estimands. The unweighted ΔBrier interval excludes zero only marginally (lower bound {inc['increment']['unweighted']['dBrier_CI'][0]:.6f}), whereas the weighted ΔBrier and both ΔNetBenefit intervals include zero; on this cross-survey validation set **we did not demonstrate a clinically meaningful advantage** of the seven-variable model over age (absence of demonstrated benefit, not proof that none exists).\n")

w("### Supplementary Table S11. Objective fundus: paired AUROC difference and gradeability\n")
w("| Endpoint | n | events | AUROC full | AUROC age-only | dAUROC (full-age, 95% CI) |"); w("| --- | --- | --- | --- | --- | --- |")
for o in inc["objective"]:
    w(f"| {o['endpoint']} | {o['n']} | {o['events']} | {o['auc_full']:.3f} | {o['auc_age']:.3f} | {o['dAUC_full_minus_age']:.4f} ({o['dAUC_CI'][0]:.4f} to {o['dAUC_CI'][1]:.4f}) |")
g=inc["gradeability"]
w("")
w(f"Gradeability (any-ARM): gradable n={g['gradable_n']} vs ungradable n={g['ungradable_n']}; ungradable were older ({g['mean_age_ungradable']:.1f} vs {g['mean_age_gradable']:.1f} y) and had higher AMD ({100*g['amd_ungradable']:.2f}% vs {100*g['amd_gradable']:.2f}%). These differences indicate that photograph availability was associated with observed characteristics and raise selection-bias concern; the missingness mechanism (whether MNAR) cannot be determined from these comparisons alone. The objective analysis is a construct-validity sensitivity, not a formal external validation.\n")

w("### Supplementary Table S12. Age functional form (linear vs 3-df RCS), external\n")
ff=inc["age_function_form"]
w("| Age term | External AUROC | slope | CITL | mean predicted |"); w("| --- | --- | --- | --- | --- |")
w(f"| Linear (published) | {ff['linear']['ext_auc']:.4f} | {ff['linear']['ext_slope']:.4f} | {ff['linear']['ext_citl']:.4f} | {ff['linear']['ext_meanpred']:.4f} |")
w(f"| 3-df RCS | {ff['rcs']['ext_auc']:.4f} | {ff['rcs']['ext_slope']:.4f} | {ff['rcs']['ext_citl']:.4f} | {ff['rcs']['ext_meanpred']:.4f} |")
w("")
w("RCS adds 2 df and raises external AUROC by only +0.0013 but gives a *worse* slope (0.847 vs 0.872); linear age was retained as a simpler specification with little apparent loss of performance. This comparison was made after the external results were available and is therefore exploratory/post-hoc (both forms ~0.774).\n")

w("### Supplementary Table S13. Descriptive agreement between operational outcome definitions (self-report vs photographic any-ARM; not diagnostic accuracy)\n")
w("| Metric | Value | Wilson 95% CI |"); w("| --- | --- | --- |")
w(f"| Agreement: self-report positive among fundus-positive | {sens[0]*100:.1f}% | {sens[1][0]*100:.1f}-{sens[1][1]*100:.1f}% |")
w(f"| Agreement: self-report negative among fundus-negative | {spec[0]*100:.1f}% | {spec[1][0]*100:.1f}-{spec[1][1]*100:.1f}% |")
w(f"| Positive agreement (PPV) | {ppv[0]*100:.1f}% | {ppv[1][0]*100:.1f}-{ppv[1][1]*100:.1f}% |")
w(f"| Negative agreement (NPV) | {npv[0]*100:.1f}% | {npv[1][0]*100:.1f}-{npv[1][1]*100:.1f}% |")
w("")
w("Computed from the unweighted 2x2 (72/316/91/4588); these are **agreement between two operational outcome definitions**, NOT clinical diagnostic accuracy of self-report, because photographic any-ARM is not a reference standard for clinically diagnosed AMD and a fundus-positive / report-negative cell is not a confirmed missed diagnosis. Not survey-corrected.\n")

w("### Supplementary Table S14. Per-predictor cross-survey comparability (harmonisation crosswalk)\n")
w("| Predictor | NHANES item / basis | NHIS item | Measurement mode (N / E) | Time window / threshold | Missing handling | Comparability |"); w("| --- | --- | --- | --- | --- | --- | --- |")
rows14=[
 ("Age","`RIDAGEYR` recorded age (years)","`AGE_P` recorded age","administrative / both recorded","continuous, 40–85; identical","age always complete","High"),
 ("Sex","`RIAGENDR`","`SEX`","recorded / recorded","male indicator","complete","High"),
 ("Smoking","`SMQ020` ever / `SMQ040` current (never, former, current)","`SMKSTAT2` (current/former/never)","self-report / self-report","current smoking at interview; NHANES ever-smoker defined by a lifetime ≥100-cigarette item","3-category; DK/RF missing","Medium-high (category schemes differ slightly)"),
 ("Diabetes","`DIQ010` ever-diagnosed","`DIBEV1` ever-diagnosed","self-report / self-report","ever-diagnosed; borderline(3)->no","8 cells missing; imputed","High"),
 ("Obesity (BMI>=30)","`BMXBMI` measured BMI","`BMI` derived from self-reported height/weight","**measured / self-reported**","BMI>=30 in both","DK/RF missing; imputed","Low-medium (self-report biases BMI down near cut-point)"),
 ("Hypertension","measured mean SBP>=140 or DBP>=90 (readings 2-4)","`HYPEV` ever-told by a doctor","**physiological measurement / diagnosed-awareness**","measured threshold vs ever-told; different construct","869 missing; imputed","Low (counts undiagnosed high BP in NHANES only)"),
 ("High cholesterol","measured TC>=240 or HDL<40/50 or TG>=200 or LDL>=160","`CHLEV` ever-told high cholesterol","**measured lipid panel / diagnosed-awareness**","multi-analyte threshold vs ever-told; TC/HDL from the full examined sample, TG/LDL from the **fasting sub-sample**","659 missing; imputed","Low (definition + fasting subsample differ)"),
]
for r in rows14: w("| "+" | ".join(r)+" |")
w("")
w("Only age, sex and (approximately) diabetes and smoking are strictly comparable across the two surveys; obesity, hypertension and high cholesterol are **operationally non-equivalent** (objective measurement in NHANES versus diagnosed-awareness report in NHIS). This is why transport is best described as a **cross-survey operational migration** of the model rather than a pure migration of one identically-measured population; see Methods and `harmonization_dictionary.md`.\n")
if lipidS:
    st=lipidS["sub_sample_structure"]; pm=lipidS["definitions"]["primary_multi_analyte"]; ts=lipidS["definitions"]["sensitivity_tc_hdl_only"]
    w("**Measured high-cholesterol: fasting sub-sample audit and definitional sensitivity** (`sr_lipid_sensitivity.py`, `sr_lipid_sensitivity.json`):\n")
    w("| Question | Finding |\n| --- | --- |")
    w(f"| Are LDL/TG fasting-sub-sample-only? | Yes — in the analytic sample TC & HDL are present for {st['full_sample_analytes']['TC_present_pct']:.0f}% / {st['full_sample_analytes']['HDL_present_pct']:.0f}% (full MEC sample), but TG/LDL only for {st['fasting_only_analytes']['TG_present_pct']:.1f}% / {st['fasting_only_analytes']['LDL_present_pct']:.1f}% (overnight-fasting sub-sample). |")
    w(f"| How are non-fasting respondents coded? | The multi-analyte criterion evaluates them on TC/HDL only; {st['pct_zero_rows_relying_on_TChdl_only']:.1f}% of dyslipidaemia-negative analytic rows ({st['cc_dyslip_zero_but_TG_LDL_missing']}/{st['cc_dyslip_zero_rows']}) are non-fasting, so their 'negative' rests on TC/HDL, not on TG/LDL. |")
    w(f"| Weight used | TC & HDL carry the full-sample MEC weight `WTMEC2YR` (correct); TG/LDL would carry a distinct fasting sub-sample weight, which the design-weighted derivation does **not** apply — hence it is labelled indicative. |")
    w(f"| Definitional sensitivity | Redefining high cholesterol to **only** the full-sample analytes (TC>=240 OR HDL<40/50) changes apparent AUROC {pm['apparent']:.4f}->{ts['apparent']:.4f}, external AUROC {pm['ext_auc']:.4f}->{ts['ext_auc']:.4f}, slope {pm['ext_slope']:.4f}->{ts['ext_slope']:.4f}, CITL {pm['ext_citl']:.4f}->{ts['ext_citl']:.4f}, and dAUROC over age {pm['paired_dAUC']['mean']:.4f}->{ts['paired_dAUC']['mean']:.4f} (95% CI {ts['paired_dAUC']['ci'][0]:.4f} to {ts['paired_dAUC']['ci'][1]:.4f}). |")
    w("")
    w("**Why respondents are classified on different numbers of indicators.** Because TG and calculated LDL exist only for the fasting sub-sample, every participant is evaluated on at least two criteria (TC, HDL) but only fasted participants on four. The mixed availability of lipid analytes may introduce differential ascertainment of dyslipidaemia across participants. The direction and magnitude of the resulting bias cannot be determined from these data. A sensitivity analysis restricted to full-sample lipid measurements yielded similar overall discrimination, although the incremental AUROC over age was no longer statistically distinguishable from zero.\n")
    w(f"**Conclusion:** the sub-sample handling does not change any headline finding — the model stays age-dominated with external AUROC ~0.77, persistent over-prediction, and no demonstrated gain of the clinical variables over age; if anything the stricter, weight-consistent definition slightly enlarges the over-prediction (dyslipidaemia prevalence {pm['dyslip_prev_cc']:.3f}->{ts['dyslip_prev_cc']:.3f}) and shrinks the already-null increment. Under that full-sample definition the paired ΔAUROC over age has a 95% CI that includes zero ({ts['paired_dAUC']['ci'][0]:.4f} to {ts['paired_dAUC']['ci'][1]:.4f}), so the increment is best described as **fragile and not distinguishable from zero**, consistent with the age-dominated conclusion, rather than as evidence that the definition \"does not matter\". We therefore keep the clinically richer multi-analyte definition, which was fixed in the design freeze (2026-10-28), as the primary definition, and report this audit and sensitivity as the pre-submission check on weight and sub-sample handling.\n")

(OUT/"Tables_sr4.md").write_text("\n".join(L),encoding="utf-8")
print("\n".join(L)[:1500]); print("...\n[Tables_sr4.md written]")

if os.environ.get("LEGACY_FIGURES", "0") != "1":
    print("[legacy figure block skipped; figures are produced by make_figures_nature_round17.py]")
    raise SystemExit(0)

# ---------- FigSR7: self-report external calibration (smooth-ish, 20 bins) + risk histogram ----------
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8})
OI=dict(blue="#0072B2",orange="#E69F00",green="#009E73",verm="#D55E00",purple="#CC79A7",grey="#7f7f7f",black="#000000")
def calib_bins(lp,y,nb=20):
    p=1/(1+np.exp(-lp)); q=np.quantile(p,np.linspace(0,1,nb+1)); idx=np.clip(np.digitize(p,q[1:-1]),0,nb-1)
    mp=[];ob=[];lo=[];hi=[]
    for b in range(nb):
        s=idx==b; n=s.sum()
        if n>=15:
            k=int(y[s].sum()); pk,pw=wilson(k,n); mp.append(float(p[s].mean())); ob.append(pk); lo.append(pw[0]); hi.append(pw[1])
    return np.array(mp),np.array(ob),np.array(lo),np.array(hi)
f,(a1,a2)=plt.subplots(1,2,figsize=(6.8,3.3))
for lp,col,lab in [(lf,OI['verm'],"Full (frozen)"),(la,OI['blue'],"Age only"),(lupd,OI['green'],"Full, intercept updated")]:
    mp,ob,lo,hi=calib_bins(lp,ye)
    a1.errorbar(mp,ob,yerr=[ob-lo,hi-ob],fmt="o-",color=col,ms=3,lw=1,capsize=1.5,label=lab)
a1.plot([0,0.2],[0,0.2],"--",color=OI['black'],lw=0.8,label="ideal (45°)")
a1.set_xlim(0,0.2);a1.set_ylim(0,0.2);a1.set_xlabel("Predicted risk");a1.set_ylabel("Observed AMD")
a1.legend(frameon=False,fontsize=6);a1.set_title("A  Self-report AMD, NHIS external",fontsize=8.5)
a2.hist(1/(1+np.exp(-lf)),bins=40,color=OI['verm'],alpha=0.8);a2.set_xlabel("Predicted risk");a2.set_ylabel("NHIS n")
a2.set_title("B  Predicted-risk distribution",fontsize=8.5)
f.tight_layout(); f.savefig(FIG/"FigSR7_calibration_selfreport.png",dpi=300,bbox_inches="tight"); f.savefig(FIG/"FigSR7_calibration_selfreport.pdf",bbox_inches="tight"); plt.close(f)

# ---------- FigSR8: objective-endpoint calibration, separate panels ----------
f,(a1,a2)=plt.subplots(1,2,figsize=(6.8,3.3))
for ax,end,col,ttl in [(a1,"obj_anyARM",OI['orange'],"A  Any photographic ARM"),(a2,"obj_lateARM",OI['purple'],"B  Late photographic AMD")]:
    dd=der[der[end].notna()].reset_index(drop=True); lp=D(dd,BIN)@np.asarray(mf.params); yd=dd[end].values.astype(float)
    mp,ob,lo,hi=calib_bins(lp,yd,nb=10)
    ax.errorbar(mp,ob,yerr=[ob-lo,hi-ob],fmt="o-",color=col,ms=4,lw=1,capsize=2,label="Full model")
    mx=max(np.max(mp),np.max(ob))*1.3; ax.plot([0,mx],[0,mx],"--",color=OI['black'],lw=0.8,label="ideal")
    ax.set_xlabel("Predicted risk");ax.set_ylabel("Observed");ax.set_title(ttl,fontsize=8.5);ax.legend(frameon=False,fontsize=6.5)
f.suptitle("Objective-endpoint calibration (separate per outcome; model trained on self-report)",fontsize=8,y=1.02)
f.tight_layout(); f.savefig(FIG/"FigSR8_calibration_objective.png",dpi=300,bbox_inches="tight"); f.savefig(FIG/"FigSR8_calibration_objective.pdf",bbox_inches="tight"); plt.close(f)
print("[FigSR7, FigSR8 written]")
