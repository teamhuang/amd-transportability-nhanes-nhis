## Main tables

### Table 1. Derivation vs external population (see Tables.md)

### Table 2. Nested model comparison

| Model | NHANES apparent AUROC | NHIS-pooled AUROC | NHIS Brier | NHIS calib slope | NHIS calib intercept |
| --- | --- | --- | --- | --- | --- |
| Age only (RCS) | 0.7831 | 0.7701 | 0.0332 | 0.9106 | -0.1882 |
| Age + sex | 0.7884 | 0.7695 | 0.0333 | 0.8791 | -0.2258 |
| Mode-clean (drop obesity/HTN/lipid) | 0.789 | 0.772 | 0.03332 | 0.8686 | -0.2236 |
| Full (7 predictors) | 0.7908 | 0.775 | 0.0334 | 0.8468 | -0.2373 |

Incremental discrimination of the full model over age alone on external data: dAUC = 0.0049 (95% bootstrap CI 0.0017-0.0086).

### Table 3. External transportability, frozen full model

| External set | n | AMD events | Prevalence | AUROC (95% CI)* | Calib slope | Calib intercept | Mean pred | Net benefit @prev |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NHIS_pooled | 38683 | 1392 | 0.0360 | 0.775 (0.762-0.787) | 0.848 | -0.233 | 0.0442 | 0.0148 |
| NHIS_2016 | 21364 | 722 | 0.0338 | 0.775 (0.758-0.793) | 0.845 | -0.296 | 0.0440 | 0.0140 |
| NHIS_2017 | 17319 | 670 | 0.0387 | 0.774 (0.755-0.791) | 0.852 | -0.160 | 0.0446 | 0.0161 |
| NHIS_region_Northeast | 6705 | 226 | 0.0337 | 0.783 (0.756-0.812) | 0.876 | -0.354 | 0.0462 | 0.0142 |
| NHIS_region_Midwest | 8754 | 351 | 0.0401 | 0.786 (0.762-0.809) | 0.890 | -0.161 | 0.0462 | 0.0173 |
| NHIS_region_South | 13941 | 467 | 0.0335 | 0.761 (0.739-0.783) | 0.798 | -0.287 | 0.0433 | 0.0127 |
| NHIS_region_West | 9283 | 348 | 0.0375 | 0.777 (0.753-0.803) | 0.859 | -0.141 | 0.0425 | 0.0159 |


*95% CI from individual bootstrap; PSU-cluster-bootstrap CI for the pooled set: AUROC [0.762, 0.7884], slope [0.7964, 0.8999], intercept [-0.2979, -0.1841].

### Table 4. Recalibration tiers, cross-year recalibration, weighting

| Scenario (NHIS pooled unless noted) | AUROC | Brier | Mean predicted | Observed |
| --- | --- | --- | --- | --- |
| Frozen, as-is | 0.7748 | 0.03337 | 0.0443 | 0.036 |
| Intercept-only recalibration | 0.7748 | 0.03306 | 0.036 | 0.036 |
| Intercept+slope recalibration | 0.7748 | 0.033 | 0.036 | 0.036 |
| 2017 frozen (no update) | - | 0.03548 | 0.0447 | 0.0387 |
| 2017 recalibrated with 2016 intercept | - | 0.03525 | 0.0342 | 0.0387 |
| Design-weighted derivation (external AUROC) | 0.7756 | - | - | - |


### Table 5. Objective fundus-grading sensitivity

| Objective endpoint (NHANES) | n | Events | Prevalence | AUROC full | AUROC age-only | Full calib slope | Full mean pred vs obs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| any AMD-related maculopathy | 5067 | 388 | 0.07657 | 0.7483 | 0.7628 | 0.8649 | 0.034 vs 0.077 |
| late AMD | 5100 | 47 | 0.00922 | 0.8937 | 0.901 | 1.6571 | 0.034 vs 0.009 |


| Cross-tab (n=5067) | Count |
| --- | --- |
| Fundus-positive AND self-report-positive | 72 |
| Fundus-positive, self-report-negative | 316 |
| Fundus-negative, self-report-positive | 91 |
| Both negative | 4588 |


Gradeability: 5560 of age>=40 had valid worst-eye grading; 12.2% of complete cases lacked a gradeable fundus photograph.

## Supplementary tables

### Supplementary Table S1. Frozen risk equation coefficients (logit)

| Term | Coefficient (logit) |
| --- | --- |
| Intercept | -8.7706 |
| age (linear RCS term) | 0.08182 |
| R1(age) nonlinear | -4e-05 |
| R2(age) nonlinear | 0.00015 |
| male | -0.46834 |
| smoke_former | 0.33287 |
| smoke_current | 0.14095 |
| obesity | 0.20144 |
| hypertension | -0.01854 |
| diabetes | 0.11875 |
| dyslipidemia | 0.27044 |

Age RCS: Harrell 3 df, knots=[41,53,65,80] . Reference = male=0 (female), non-smoker, non-obese, normotensive, non-diabetic, normal lipids. Worked examples — 70y male (former smoker, obese, HTN, high-chol, non-diab): lp=-2.9847, p=0.0481; 55y healthy female: p=0.0124.

### Supplementary Table S2. Missingness

| Variable / metric | Value |
| --- | --- |
| Missing obesity | 6.1% |
| Missing hypertension | 12.4% |
| Missing dyslipidemia | 9.4% |
| Missing diabetes/smoking | ~0.1% |
| Any predictor missing (age>=40) | 17.7% |
| AMD prevalence, included (complete-case) | 3.71% |
| AMD prevalence, all age>=40 | 4.12% |
| NHIS any predictor missing | 4.3% |


### Supplementary Table S3. Subgroup performance (NHIS)

| NHIS subgroup | n? | AUROC | Calib slope | Calib intercept | Mean pred vs obs |
| --- | --- | --- | --- | --- | --- |
| age 40-54 |  | 0.6504 | 1.136 | 0.0229 | 0.008 vs 0.008 |
| age 55-64 |  | 0.5779 | 0.7304 | 0.0022 | 0.019 vs 0.019 |
| age 65-74 |  | 0.5675 | 0.5355 | -0.1495 | 0.045 vs 0.039 |
| age 75+ |  | 0.6326 | 0.9037 | -0.3773 | 0.154 vs 0.112 |
| sex male |  | 0.7651 | 0.8313 | -0.0804 | 0.032 vs 0.030 |
| sex female |  | 0.7794 | 0.8836 | -0.3272 | 0.055 vs 0.041 |


### Supplementary Table S4. Decision-curve net benefit, prespecified thresholds

| Threshold | Model net benefit | Treat-all net benefit | Flagged /1000 | Cases caught /1000 | Sensitivity | Cases missed /1000 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.02 | 0.02102 | 0.01631 | 501 | 30 | 0.851 | 5 |
| 0.05 | 0.01091 | -0.01475 | 243 | 22 | 0.626 | 13 |
| 0.08 | 0.00607 | -0.04784 | 162 | 18 | 0.515 | 17 |
| 0.1 | 0.00412 | -0.07113 | 131 | 16 | 0.468 | 19 |


### Supplementary Table S5. Age-only reference risk curve

| Age (y) | Predicted AMD risk (reference profile) |
| --- | --- |
| 45 | 0.0061 |
| 50 | 0.0089 |
| 55 | 0.0124 |
| 60 | 0.0166 |
| 65 | 0.0231 |
| 70 | 0.0355 |
| 75 | 0.0587 |
| 80 | 0.0993 |
| 85 | 0.164 |

