### Table 1. Study population characteristics

| Characteristic | NHANES (derivation) | NHIS 2016-2017 (external) |
| --- | --- | --- |
| n (>=40, complete-case) | 5773 | 38683 |
| Age, mean (SD) | 59.7 (12.6) | 61.2 (12.5) |
| Female (%) | 49.8 | 54.1 |
| Former smoker (%) | 31.9 | 30.3 |
| Current smoker (%) | 20.0 | 15.3 |
| Obesity BMI>=30 (%) | 37.4 | 32.5 |
| Hypertension (%) | 25.7 | 47.1 |
| Diabetes (%) | 15.8 | 14.5 |
| High cholesterol (%) | 44.6 | 42.2 |
| Self-reported AMD (%) | 3.71 | 3.60 |

### Table 2. NHANES derivation model (risk equation) — adjusted effects

| Predictor | Adjusted OR (95% CI) | p (raw) | p (Holm) |
| --- | --- | --- | --- |
| Male (ref: female) | 0.63 (0.46-0.85) | 0.00229 | 0.0229 |
| Former smoker (ref: never) | 1.39 (1.01-1.90) | 0.0402 | 0.362 |
| Current smoker (ref: never) | 1.15 (0.71-1.87) | 0.574 | 1 |
| Obesity BMI>=30 | 1.22 (0.90-1.66) | 0.207 | 1 |
| Hypertension (measured >=140/90) | 0.98 (0.73-1.33) | 0.908 | 1 |
| Diabetes (ever dx) | 1.12 (0.79-1.60) | 0.521 | 1 |
| High cholesterol / dyslipidemia | 1.31 (0.98-1.75) | 0.0642 | 0.513 |
| Age (RCS spline); OR 80 vs 55 y | 11.3 (vs 55) | - | dominant |

### Table 3. External validation (transportability) performance by set

| External set | n | AMD events | Prevalence | AUROC (95% CI) | Calib slope | Calib intercept | Brier | Mean pred prob | Net benefit @prev |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NHIS_pooled | 38683 | 1392 | 0.0360 | 0.775 (0.762-0.787) | 0.848 | -0.233 | 0.0334 | 0.0442 | 0.0148 |
| NHIS_2016 | 21364 | 722 | 0.0338 | 0.775 (0.758-0.793) | 0.845 | -0.296 | 0.0317 | 0.0440 | 0.0140 |
| NHIS_2017 | 17319 | 670 | 0.0387 | 0.774 (0.755-0.791) | 0.852 | -0.160 | 0.0355 | 0.0446 | 0.0161 |
| NHIS_region_Northeast | 6705 | 226 | 0.0337 | 0.783 (0.756-0.812) | 0.876 | -0.354 | 0.0315 | 0.0462 | 0.0142 |
| NHIS_region_Midwest | 8754 | 351 | 0.0401 | 0.786 (0.762-0.809) | 0.890 | -0.161 | 0.0365 | 0.0462 | 0.0173 |
| NHIS_region_South | 13941 | 467 | 0.0335 | 0.761 (0.739-0.783) | 0.798 | -0.287 | 0.0317 | 0.0433 | 0.0127 |
| NHIS_region_West | 9283 | 348 | 0.0375 | 0.777 (0.753-0.803) | 0.859 | -0.141 | 0.0343 | 0.0425 | 0.0159 |

### Table 4. Recalibration tiers and weighted sensitivity (NHIS pooled)

| Scenario | AUROC | Brier | Mean predicted prob |
| --- | --- | --- | --- |
| Frozen model, as-is | 0.7748 | 0.03337 | 0.0443 |
| Intercept-only recalibration | 0.7748 | 0.03306 | 0.036 |
| Intercept+slope recalibration | 0.7748 | 0.033 | 0.036 |
| Design-weighted derivation (sensitivity) | 0.7756 | - | - |

### Table 5. Objective fundus-grading sensitivity (NHANES)

| Endpoint (NHANES only) | n | Events | Prevalence | AUROC of frozen self-report score (95% CI) |
| --- | --- | --- | --- | --- |
| Self-reported AMD (derivation) | 5773 | 214 | 0.0371 | 0.782 (bootstrap-corrected) |
| Objective any-ARM (fundus) | 5067 | 388 | 0.0766 | 0.749 (0.722-0.772) |
| Objective late-ARM (fundus) | 5100 | 47 | 0.00922 | 0.894 (0.851-0.933) |


| Concordance self-report vs fundus any-ARM | Count |
| --- | --- |
| Both positive | 72 |
| Fundus-positive only (undiagnosed) | 316 |
| Self-report-positive only | 91 |
| Neither | 4588 |
| Total (both measured) | 5067 |

