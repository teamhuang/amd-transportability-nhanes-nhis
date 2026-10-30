## Round-4 supplementary tables

### Supplementary Table S7. Endpoint definitions verified against the live CDC codebook and the raw NCHS transport files (audit log: `endpoint_code_audit.md`)

| Endpoint | Survey / variable | Verified wording / coding | Our operationalisation |
| --- | --- | --- | --- |
| Self-reported AMD (outcome) | NHANES `VIQ310` (2005–06, 2007–08) | "…ever been told by an **eye doctor** that … age-related macular degeneration?"; both sexes 40–150 y; 1=Yes, 2=No, 7=Refused (0 occurrences), 9=Don't know | `VIQ310==1`; 7/9 missing |
| Self-reported AMD (outcome) | NHIS `VIM_MDEV` (2016, 2017) | "EVER been told you had macular degeneration"; sample adults ≥18 y; 1=Yes, 2=No, 7=Refused, 8=Not ascertained (0 occurrences), 9=Don't know | `VIM_MDEV==1`; 7/9 missing |
| Photographic any-ARM | NHANES `OPDUARM` (2005–06, 2007–08) | "ARM, 3 severity levels, worse eye"; both sexes 40–150 y; 0=No ARM, 1=Early ARM (drusen and/or pigmentary abnormalities), 2=Late ARM (exudative ARM and/or geographic atrophy) | `OPDUARM>=1` (early+late photographic ARM) |
| Photographic late AMD | NHANES `OPDULARM` (2005–06, 2007–08) | "Late ARM, worse eye"; 0=No, 1=Yes (exudative ARM and/or geographic atrophy) | `OPDULARM==1` |

**Verification status (completed before submission, October 2026):** every label and code was read both from the raw `.xpt`/`.csv` transports and from the live CDC/NCHS codebook pages; counts agree **cell-by-cell and exactly** in both NHANES cycles and both NHIS years (`analysis/scripts/audit_codebook_round12.py`; machine-readable record `analysis/output/sr_codebook_audit.json`). **Caveat:** CDC publishes self-reported AMD (`VIQ310`, `VIM_MDEV`) and examination-based AMD (`OPDUARM`) as **separate case definitions** and makes **no equivalence statement** between them, and a third claims-based "diagnosed AMD" definition uses ICD-9/ICD-10 codes in claims or EHR data (CDC Vision Health Surveillance case definitions, https://www.cdc.gov/vision-health-data/case-definitions/age-related-macular-degeneration-amd.html). The NHANES and NHIS self-report items also differ in wording (eye doctor vs any doctor), an outcome-ascertainment difference across surveys.

### Supplementary Table S8. Sample flow and complete-case vs multiple-imputation comparison

| Step / metric | NHANES n (or value) |
| --- | --- |
| Age >=40 with valid AMD outcome | 7011 |
| Complete predictors (analytic) | 5773 |
| Excluded, incomplete predictors | 1238 (17.7%) |
| Fundus-gradable any-ARM | 5067 |
| Fundus-gradable late ARM | 5100 |
| AMD prevalence: included vs excluded vs all | 3.71% vs 6.06% vs 4.12% |

**Complete-case vs multiple-imputation comparison (same S8):**

**The missing-data sensitivity uses chained-equations multiple imputation with variable-type-appropriate conditional models** (logistic for the four binary health indicators; multinomial logistic for the three-category smoking variable never/former/current), Bayesian parameter draws `beta* ~ N(beta_hat, Cov_hat)`, m = 20 imputations x 15 sweeps, outcome-inclusive and cycle-adjusted. Out-of-range and illegal-combination imputations are now structurally impossible (categorical sampling), so no clipping/rounding is used. (An earlier linear-Gaussian `IterativeImputer`+`BayesianRidge` version that clipped 12,530 out-of-range cells was superseded as statistically inappropriate for these variable types.)

| Quantity | Complete-case (primary) | Multiple imputation (m=20) |
| --- | --- | --- |
| External AUROC | 0.7737 | 0.7744 (range 0.7723-0.7756, between-imp SD 0.0009) |
| External calibration slope | 0.8723 | 0.9018 |
| External CITL intercept | -0.2131 | -0.2204 |
| Age coefficient (logit), Rubin-pooled | 0.0945 | 0.0932 (95% CI 0.0801 to 0.1064, total SE 0.0067) |
| Imputation validity (m x n) | n/a | 0 out-of-range cells; 0 illegal smoking combinations (both structurally prevented) |

**MI diagnostics (observed vs imputed distributions):**

| Variable | Model type | % missing | Observed category % | Imputed-overall category % |
| --- | --- | --- | --- | --- |
| obesity | logistic | 6.12 | 37.4% positive | 37.7% positive |
| hypertension | logistic | 12.39 | 25.7% positive | 26.9% positive |
| dyslipidemia | logistic | 9.40 | 44.6% positive | 44.4% positive |
| diabetes | logistic | 0.11 | 15.8% positive | 16.4% positive |
| smoking (3-category) | multinomial-logit | 0.10 | never 48.5 / former 31.8 / current 19.6 | never 48.5 / former 31.9 / current 19.6 |

MI recovers the 1,238 excluded participants and leaves the core conclusions unchanged (age-dominated, external AUROC ~0.774, persistent over-prediction). Observed and imputed category proportions are close (S8 diagnostics block), indicating the type-appropriate conditional models reproduce the empirical distribution rather than distorting it. Complete-case selection is **not** negligible (excluded AMD 6.06% vs included 3.71%) but does not overturn the findings. The between-imputation **range** (0.7723-0.7756) reflects imputation variability; formal inference for a coefficient uses Rubin's rules (total SE, Barnard-Rubin degrees of freedom), and the Barnard-Rubin correction here is negligible (very large residual df) so the reported age-coefficient interval is effectively a normal interval.

**Independent code re-check of the MICE (empirical convergence across seeds x sweeps):** 
re-running the sampler over sweep counts 5/15/30 x three random seeds at m=6 gives a pooled external-AUROC grid mean **0.7746** with total spread **0.0006** across the whole grid (published single-run value 0.7744); across 4,500 conditional fits there were **0** binary and **0** smoking fallbacks to the marginal (no singular-Hessian / perfect-separation event occurred), and imputed values were valid 0/1 or a single smoking level in every draw. Run log archived as `sr_mi_stability_runlog.txt`. The absence of a formal MICE trace is why we report stability across cycles+seeds rather than claiming a convergence statistic.

### Supplementary Table S9. Complex-survey external inference (intervals are method-specific, not one common design-based estimator)

| Quantity (NHIS) | Estimate | 95% interval | Method (see caveat) |
| --- | --- | --- | --- |
| Weighted AMD prevalence (pooled) | 2.88% | 2.67-3.09% | Taylor linearisation (strata+PSU) |
| Weighted AMD prevalence (2016) | 2.68% | 2.40-2.96% | Taylor linearisation (strata+PSU) |
| Weighted AMD prevalence (2017) | 3.07% | 2.76-3.38% | Taylor linearisation (strata+PSU) |
| Weighted calibration slope (pooled) | 0.903 | 0.836-0.971 | Weighted GLM, PSU-clustered SE |
| Weighted calibration slope (pooled) | 0.903 | 0.842-0.969 | PSU bootstrap B=1000 (seed 1) |
| Weighted CITL (pooled) | -0.237 | -0.310--0.165 | Weighted GLM, PSU-clustered SE |
| Weighted AUROC (pooled) | 0.787 | 0.772-0.803 | PSU bootstrap B=1000 |

Bootstrap seed-to-seed stability (max CI-endpoint difference across seeds 20261028 and 777): 0.0066. Weighted AUROC (0.787) differs from unweighted (0.774) because it targets a different population, not because the model is better.

### Supplementary Table S10. Full-vs-age-only increment (paired bootstrap), by estimand

| Estimand | ΔAUROC (full−age) (95% CI) | ΔBrier (full−age) (95% CI) | ΔNetBenefit@5% (full−age) (95% CI) |
| --- | --- | --- | --- |
| unweighted | 0.0037 (0.0008 to 0.0068) | 0.000111 (0.000003 to 0.000217) | -0.00003 (-0.00050 to 0.00043) |
| weighted | 0.0048 (0.0015 to 0.0082) | 0.000049 (-0.000044 to 0.000138) | -0.00012 (-0.00060 to 0.00039) |

All three deltas are **full-model minus age-only**. A **positive ΔBrier means the full model has a HIGHER (worse) Brier score** than age only; a positive ΔAUROC/ΔNetBenefit favours the full model. The AUROC gain is statistically detectable but tiny in both estimands. The unweighted ΔBrier interval excludes zero only marginally (lower bound 0.000003), whereas the weighted ΔBrier and both ΔNetBenefit intervals include zero; on this cross-survey validation set **we did not demonstrate a clinically meaningful advantage** of the seven-variable model over age (absence of demonstrated benefit, not proof that none exists).

### Supplementary Table S11. Objective fundus: paired AUROC difference and gradeability

| Endpoint | n | events | AUROC full | AUROC age-only | dAUROC (full-age, 95% CI) |
| --- | --- | --- | --- | --- | --- |
| obj_anyARM | 5067 | 388 | 0.751 | 0.763 | -0.0117 (-0.0187 to -0.0052) |
| obj_lateARM | 5100 | 47 | 0.894 | 0.901 | -0.0069 (-0.0223 to 0.0080) |

Gradeability (any-ARM): gradable n=5067 vs ungradable n=706; ungradable were older (64.6 vs 59.1 y) and had higher AMD (7.22% vs 3.22%). These differences indicate that photograph availability was associated with observed characteristics and raise selection-bias concern; the missingness mechanism (whether MNAR) cannot be determined from these comparisons alone. The objective analysis is a construct-validity sensitivity, not a formal external validation.

### Supplementary Table S12. Age functional form (linear vs 3-df RCS), external

| Age term | External AUROC | slope | CITL | mean predicted |
| --- | --- | --- | --- | --- |
| Linear (published) | 0.7737 | 0.8723 | -0.2131 | 0.0436 |
| 3-df RCS | 0.7750 | 0.8468 | -0.2373 | 0.0444 |

RCS adds 2 df and raises external AUROC by only +0.0013 but gives a *worse* slope (0.847 vs 0.872); linear age was retained as a simpler specification with little apparent loss of performance. This comparison was made after the external results were available and is therefore exploratory/post-hoc (both forms ~0.774).

### Supplementary Table S13. Descriptive agreement between operational outcome definitions (self-report vs photographic any-ARM; not diagnostic accuracy)

| Metric | Value | Wilson 95% CI |
| --- | --- | --- |
| Agreement: self-report positive among fundus-positive | 18.6% | 15.0-22.7% |
| Agreement: self-report negative among fundus-negative | 98.1% | 97.6-98.4% |
| Positive agreement (PPV) | 44.2% | 36.8-51.8% |
| Negative agreement (NPV) | 93.6% | 92.8-94.2% |

Computed from the unweighted 2x2 (72/316/91/4588); these are **agreement between two operational outcome definitions**, NOT clinical diagnostic accuracy of self-report, because photographic any-ARM is not a reference standard for clinically diagnosed AMD and a fundus-positive / report-negative cell is not a confirmed missed diagnosis. Not survey-corrected.

### Supplementary Table S14. Per-predictor cross-survey comparability (harmonisation crosswalk)

| Predictor | NHANES item / basis | NHIS item | Measurement mode (N / E) | Time window / threshold | Missing handling | Comparability |
| --- | --- | --- | --- | --- | --- | --- |
| Age | `RIDAGEYR` recorded age (years) | `AGE_P` recorded age | administrative / both recorded | continuous, 40–85; identical | age always complete | High |
| Sex | `RIAGENDR` | `SEX` | recorded / recorded | male indicator | complete | High |
| Smoking | `SMQ020` ever / `SMQ040` current (never, former, current) | `SMKSTAT2` (current/former/never) | self-report / self-report | current smoking at interview; NHANES ever-smoker defined by a lifetime ≥100-cigarette item | 3-category; DK/RF missing | Medium-high (category schemes differ slightly) |
| Diabetes | `DIQ010` ever-diagnosed | `DIBEV1` ever-diagnosed | self-report / self-report | ever-diagnosed; borderline(3)->no | 8 cells missing; imputed | High |
| Obesity (BMI>=30) | `BMXBMI` measured BMI | `BMI` derived from self-reported height/weight | **measured / self-reported** | BMI>=30 in both | DK/RF missing; imputed | Low-medium (self-report biases BMI down near cut-point) |
| Hypertension | measured mean SBP>=140 or DBP>=90 (readings 2-4) | `HYPEV` ever-told by a doctor | **physiological measurement / diagnosed-awareness** | measured threshold vs ever-told; different construct | 869 missing; imputed | Low (counts undiagnosed high BP in NHANES only) |
| High cholesterol | measured TC>=240 or HDL<40/50 or TG>=200 or LDL>=160 | `CHLEV` ever-told high cholesterol | **measured lipid panel / diagnosed-awareness** | multi-analyte threshold vs ever-told; TC/HDL from the full examined sample, TG/LDL from the **fasting sub-sample** | 659 missing; imputed | Low (definition + fasting subsample differ) |

Only age, sex and (approximately) diabetes and smoking are strictly comparable across the two surveys; obesity, hypertension and high cholesterol are **operationally non-equivalent** (objective measurement in NHANES versus diagnosed-awareness report in NHIS). This is why transport is best described as a **cross-survey operational migration** of the model rather than a pure migration of one identically-measured population; see Methods and `harmonization_dictionary.md`.

**Measured high-cholesterol: fasting sub-sample audit and definitional sensitivity** (`sr_lipid_sensitivity.py`, `sr_lipid_sensitivity.json`):

| Question | Finding |
| --- | --- |
| Are LDL/TG fasting-sub-sample-only? | Yes — in the analytic sample TC & HDL are present for 100% / 100% (full MEC sample), but TG/LDL only for 48.5% / 47.1% (overnight-fasting sub-sample). |
| How are non-fasting respondents coded? | The multi-analyte criterion evaluates them on TC/HDL only; 50.9% of dyslipidaemia-negative analytic rows (1626/3197) are non-fasting, so their 'negative' rests on TC/HDL, not on TG/LDL. |
| Weight used | TC & HDL carry the full-sample MEC weight `WTMEC2YR` (correct); TG/LDL would carry a distinct fasting sub-sample weight, which the design-weighted derivation does **not** apply — hence it is labelled indicative. |
| Definitional sensitivity | Redefining high cholesterol to **only** the full-sample analytes (TC>=240 OR HDL<40/50) changes apparent AUROC 0.7911->0.7913, external AUROC 0.7737->0.7732, slope 0.8723->0.8597, CITL -0.2131->-0.2740, and dAUROC over age 0.0036->0.0031 (95% CI -0.0003 to 0.0066). |

**Why respondents are classified on different numbers of indicators.** Because TG and calculated LDL exist only for the fasting sub-sample, every participant is evaluated on at least two criteria (TC, HDL) but only fasted participants on four. The mixed availability of lipid analytes may introduce differential ascertainment of dyslipidaemia across participants. The direction and magnitude of the resulting bias cannot be determined from these data. A sensitivity analysis restricted to full-sample lipid measurements yielded similar overall discrimination, although the incremental AUROC over age was no longer statistically distinguishable from zero.

**Conclusion:** the sub-sample handling does not change any headline finding — the model stays age-dominated with external AUROC ~0.77, persistent over-prediction, and no demonstrated gain of the clinical variables over age; if anything the stricter, weight-consistent definition slightly enlarges the over-prediction (dyslipidaemia prevalence 0.446->0.290) and shrinks the already-null increment. Under that full-sample definition the paired ΔAUROC over age has a 95% CI that includes zero (-0.0003 to 0.0066), so the increment is best described as **fragile and not distinguishable from zero**, consistent with the age-dominated conclusion, rather than as evidence that the definition "does not matter". We therefore keep the clinically richer multi-analyte definition, which was fixed in the design freeze (2026-10-28), as the primary definition, and report this audit and sensitivity as the pre-submission check on weight and sub-sample handling.
