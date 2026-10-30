## Tables

### Table 1. Derivation vs external population (unweighted and design-weighted)

| Characteristic | NHANES unweighted | NHANES weighted | NHIS unweighted | NHIS weighted |
| --- | --- | --- | --- | --- |
| N (>=40, complete-case) | 5773 |  | 38683 |  |
| Age, mean (y) | 59.75 | 56.55 | 61.23 | 59.02 |
| Female (%) | 49.8 | 52.2 | 54.1 | 51.6 |
| Current smoker (%) | 19.97 | 20.37 | 15.26 | 14.4 |
| Obesity BMI>=30 (%) | 37.4 | 36.21 | 32.53 | 32.86 |
| Hypertension (%) | 25.71 | 21.65 | 47.14 | 44.06 |
| Diabetes (%) | 15.76 | 11.37 | 14.53 | 13.81 |
| High cholesterol (%) | 44.62 | 45.09 | 42.15 | 40.47 |
| Self-reported AMD (%) | 3.71 | 3.11 | 3.6 | 2.88 |


### Table 2. Nested models (frozen, linear-age canonical)

| Model | Derivation AUROC | NHIS-pooled AUROC | NHIS slope | NHIS intercept | NHIS Brier |
| --- | --- | --- | --- | --- | --- |
| Age only | 0.7831 | 0.7701 | 0.94 | -0.1672 | 0.03314 |
| Age+sex | 0.7874 | 0.7695 | 0.903 | -0.2048 | 0.0332 |
| Mode-clean | 0.7886 | 0.7716 | 0.8947 | -0.1998 | 0.0332 |
| Full (7) | 0.7911 | 0.7737 | 0.8723 | -0.2131 | 0.03325 |


Incremental performance of the full model over age only (external; single frozen 1,000-replicate paired bootstrap, seed 20261028). Unweighted analytic-sample estimand: dAUROC 0.0037 (95% CI 0.0008-0.0068); dBrier 0.000111 (95% CI 0.000003 to 0.000217); dNetBenefit@5% -0.00003 (95% CI -0.00050 to 0.00043). Design-weighted estimand: dAUROC 0.0048 (95% CI 0.0015-0.0082); dBrier 0.000049 (95% CI -0.000044 to 0.000138). All dAUROC/dBrier/dNB values in this manuscript come from this one run; the PSU-cluster interval for the AUROC point estimate (not the difference) is given in Table 3.

### Table 3. External transportability, frozen full model

| External set (frozen full model) | n | Events | Prevalence (%) | AUROC (indiv. boot 95% CI) | Calib slope | Calib intercept | Brier | Mean pred vs obs (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pooled | 38683 | 1392 | 3.60 | 0.774 (0.761-0.786) | 0.8723 | -0.2131 | 0.03325 | 4.36 vs 3.60 |
| 2016 | 21364 | 722 | 3.38 | 0.775 (0.757-0.792) | 0.8717 | -0.2757 | 0.03152 | 4.34 vs 3.38 |
| 2017 | 17319 | 670 | 3.87 | 0.773 (0.754-0.790) | 0.8728 | -0.1411 | 0.03538 | 4.39 vs 3.87 |
| region_Northeast | 6705 | 226 | 3.37 | 0.782 (0.755-0.810) | 0.915 | -0.3307 | 0.03136 | 4.54 vs 3.37 |
| region_Midwest | 8754 | 351 | 4.01 | 0.784 (0.761-0.808) | 0.9177 | -0.1365 | 0.03639 | 4.53 vs 4.01 |
| region_South | 13941 | 467 | 3.35 | 0.760 (0.738-0.783) | 0.8183 | -0.2705 | 0.03151 | 4.28 vs 3.35 |
| region_West | 9283 | 348 | 3.75 | 0.776 (0.752-0.801) | 0.8766 | -0.1236 | 0.03427 | 4.19 vs 3.75 |


Pooled **PSU-cluster-bootstrap** 95% intervals: AUROC [0.7602, 0.7875], slope [0.8157, 0.9353], intercept [-0.2726, -0.1594] (composite-cluster bootstrap intervals — an approximation to design-based inference, not a Taylor/`svyglm` variance; the AUROC CI shown in the table is the individual-level bootstrap).

### Table 4. Recalibration / cross-year update / weighting

| Scenario | Brier | Mean predicted (%) | Observed (%) |
| --- | --- | --- | --- |
| NHIS 2017 frozen | 0.03538 | 4.39 | 3.87 |
| NHIS 2017, intercept updated from 2016 | 0.03528 | 3.42 | 3.87 |
| NHIS 2016 frozen | 0.03152 | 4.34 | 3.38 |
| NHIS 2016, intercept updated from 2017 | 0.03131 | 3.82 | 3.38 |


Frozen rows reuse the same transport-path output as Table 3, so the two tables agree exactly; intercept-updated rows use the temporal-update path. All risks are on the percentage scale used in the text. Intercept corrections: 2016 = -0.2757, 2017 = -0.1411; frozen-model AUROC 2016=0.7746, 2017=0.7726 (unchanged by intercept shift). Design-weighted derivation external AUROC = 0.7756. Note: the two cross-year directions behave differently — updating 2017 from 2016 overshoots into slight UNDER-prediction (updated 3.42% vs observed 3.87%), while updating 2016 from 2017 still OVER-predicts (updated 3.82% vs observed 3.38%); an intercept-only update zeroes calibration-in-the-large **within the sample used to estimate it** but never corrects the slope, so applied to an independent year it does not reliably restore calibration.

### Table 5. Objective fundus-grading sensitivity

| Objective endpoint (NHANES) | n | Events | Prevalence (%) | AUROC full | AUROC age-only | Full calib slope | Full mean pred vs obs (%) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| any-ARM | 5067 | 388 | 7.66 | 0.751 | 0.7628 | 0.834 | 3.42 vs 7.66 |
| late-AMD | 5100 | 47 | 0.92 | 0.8941 (CI [0.849, 0.933]) | 0.901 | 1.801 | 3.44 vs 0.92 |


| Cross-tab (n=5067) | Count |
| --- | --- |
| Fundus+ & self-report+ | 72 |
| Fundus+ , self-report- (not reporting a diagnosis) | 316 |
| Fundus- , self-report+ | 91 |
| Both - | 4588 |


Denominators: among 5773 analytic complete cases, 5067 had valid worst-eye any-ARM grading (12.23% ungradeable) and 5100 valid late-ARM grading; the two endpoints differ slightly because photograph quality/eye availability differ by grading variable. Late-AMD has only 47 events (AUROC CI [0.849, 0.933]).

## Supplementary tables

### Supplementary Table S1. Frozen risk equation (logit)

| Term | Coefficient (logit) | OR |
| --- | --- | --- |
| Intercept | -9.6317 | — |
| age | 0.0945 | 1.099 |
| male | -0.4657 | 0.628 |
| smoke_former | 0.307 | 1.359 |
| smoke_current | 0.1179 | 1.125 |
| obesity | 0.1668 | 1.182 |
| hypertension | -0.0185 | 0.982 |
| diabetes | 0.089 | 1.093 |
| dyslipidemia | 0.2686 | 1.308 |


Age entered as a LINEAR term (OR 1.0991 per year, 2.572 per decade); a 3-df restricted cubic spline was inspected only post hoc (raised external AUROC by +0.0013, 0.7750 vs 0.7737) and was **not** used to select the published linear term. Reference profile: male=0 (female), never-smoker, non-obese, normotensive, non-diabetic, normal lipids. Worked examples — {'70yM(former,obese,HTN,highchol,nonDiab)': 0.0596, '55F healthy': 0.0117}. Reproduction: with 4-dp coefficients the maximum absolute deviation in linear predictor across the cohort is 0.00223 (probability 0.000453), so `risk_calculator.py` reproduces the model. **Locations of attached artefacts:** the FULL-precision coefficients for both the linear and spline specifications and the spline knots are stored in `model_equation.json` (this table shows 4-dp values for readability); the executable calculator is `risk_calculator.py`; the per-predictor joint Wald tests with Holm-adjusted p-values are in the table below (Part B of S1).

**S1 Part B — joint Wald tests and Holm-adjusted p (Holm across the seven predictors):**

| Predictor (joint Wald) | chi-square | df | p | Holm-adjusted p |
| --- | --- | --- | --- | --- |
| Age (1 df) | 158.14 | 1 | 2.9e-36 | 0.000 |
| Male vs female (1 df) | 9.24 | 1 | 0.002 | 0.014 |
| Smoking, joint former+current (2 df) | 3.73 | 2 | 0.155 | 0.621 |
| Obesity BMI>=30 (1 df) | 1.15 | 1 | 0.284 | 0.853 |
| Hypertension (1 df) | 0.01 | 1 | 0.904 | 1.000 |
| Diabetes (1 df) | 0.25 | 1 | 0.620 | 1.000 |
| High cholesterol (1 df) | 3.37 | 1 | 0.066 | 0.332 |


### Supplementary Table S2. Missingness

| Variable / metric | NHANES |
| --- | --- |
| Missing obesity | 6.1% |
| Missing hypertension | 12.4% |
| Missing dyslipidemia | 9.4% |
| Missing diabetes or smoking | ~0.1% |
| Any predictor missing (age>=40) | 17.7% |
| AMD prevalence: included vs all age>=40 | 3.71% vs 4.12% |
| NHIS any predictor missing | 4.3% |


### Supplementary Table S3. Subgroup performance (NHIS)

| NHIS subgroup | n | Events | AUROC (95% CI) | slope | intercept | Mean pred vs obs (%) |
| --- | --- | --- | --- | --- | --- | --- |
| age 40-54 | 13037 | 110 | 0.647 (0.592-0.701) | 0.9905 | 0.2604 | 0.65 vs 0.84 |
| age 55-64 | 10166 | 198 | 0.570 (0.534-0.604) | 0.6146 | -0.0671 | 2.08 vs 1.95 |
| age 65-74 | 8898 | 345 | 0.566 (0.537-0.592) | 0.5649 | -0.3136 | 5.22 vs 3.88 |
| age 75+ | 6582 | 739 | 0.629 (0.609-0.649) | 1.0473 | -0.2633 | 14.07 vs 11.23 |
| sex male | 17737 | 532 | 0.762 (0.742-0.782) | 0.847 | -0.0656 | 3.19 vs 3.00 |
| sex female | 20946 | 860 | 0.779 (0.764-0.793) | 0.9216 | -0.2979 | 5.35 vs 4.11 |


### Supplementary Table S4. Decision-curve at the analysed thresholds (2%, 5%, 8%, 10%)

| Threshold | Flagged/1000 | Cases caught/1000 | Sensitivity | Cases missed/1000 | NB full | NB age-only | NB all | NB none |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.02 | 522 | 31 | 0.856 | 5 | 0.02077 | 0.02061 | 0.01631 | 0 |
| 0.05 | 280 | 24 | 0.667 | 12 | 0.01053 | 0.01056 | -0.01475 | 0 |
| 0.08 | 173 | 19 | 0.536 | 17 | 0.00594 | 0.0066 | -0.04784 | 0 |
| 0.1 | 136 | 17 | 0.473 | 19 | 0.00384 | 0.00422 | -0.07113 | 0 |


### Supplementary Table S5. Predicted risk by age

| Age | Full-model risk (reference profile) | Age-only model risk |
| --- | --- | --- |
| 45 | 0.0046 | 0.0055 |
| 50 | 0.0073 | 0.0086 |
| 55 | 0.0117 | 0.0136 |
| 60 | 0.0186 | 0.0213 |
| 65 | 0.0296 | 0.0334 |
| 70 | 0.0466 | 0.0519 |
| 75 | 0.0727 | 0.0798 |
| 80 | 0.1117 | 0.1208 |
| 85 | 0.1678 | 0.1788 |

### Supplementary Table S6. External calibration under two estimands (frozen full model)

| Quantity (NHIS pooled) | Unweighted (primary) | Design-weighted |
| --- | --- | --- |
| AUROC | 0.774 | 0.787 |
| Calibration slope | 0.872 | 0.903 |
| Calibration-in-the-large (logit) | -0.213 | -0.237 |
| Mean predicted risk | 4.36% | 3.57% |
| Observed prevalence | 3.60% | 2.88% |
| Over-prediction ratio (pred/obs) | 1.21 | 1.24 |
| Brier | 0.0333 | 0.0268 |

Cross-year intercept update (2016→2017): unweighted mean predicted 4.39%→3.42% against observed 3.87%; design-weighted 3.61%→2.73% against observed 3.07%. Under both estimands the update moves predicted risk closer to, but below, the observed value (Brier 0.03538→0.03528 unweighted; 0.02851→0.02840 weighted).

Decision-curve net benefit of the full model at thresholds 2/5/8/10%: unweighted 0.0208, 0.0105, 0.0059, 0.0038; design-weighted 0.0153, 0.0072, 0.0042, 0.0025 (weighted treat-all 0.0090, -0.0223, -0.0557, -0.0791). Net benefit remains above treat-all/treat-none at low thresholds under both estimands; its magnitude tracks the weighted prevalence and is therefore estimand-specific.

