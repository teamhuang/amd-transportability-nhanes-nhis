# Variable harmonization dictionary (NHANES derivation <-> NHIS external validation)
Outcome, self-reported AMD: NHANES VIQ310 (==1) <-> NHIS VIM_MDEV (==1). 7/9/77/99 -> missing.
Age: NHANES RIDAGEYR <-> NHIS AGE_P (both 40-85; NHIS not top-coded at 70 - verified: max 85, ~18% >=70).
Sex: NHANES RIAGENDR==1(male) <-> NHIS SEX==1. male indicator.
Smoking: NHANES SMQ020(ever>=100)/SMQ040(current) -> former,current (never ref) <-> NHIS SMKSTAT2 (1/2 current,3 former,4 never).
Obesity BMI>=30: NHANES BMXBMI (measured) <-> NHIS BMI (=BMI*100; 1000<=v<=6000 valid) -> self-reported. MEASUREMENT-MODE DIFFERENCE.
Hypertension: NHANES mean SBP>=140 or DBP>=90 (measured readings 2-4) <-> NHIS HYPEV==1 (ever-told). MEASUREMENT-MODE DIFFERENCE.
Diabetes: NHANES DIQ010==1 <-> NHIS DIBEV1==1 (both ever-diagnosed). 3=borderline -> not diabetic.
High cholesterol: NHANES measured TC>=240 or HDL<40/50 or TG>=200 or LDL>=160 <-> NHIS CHLEV==1 (ever-told). MEASUREMENT-MODE DIFFERENCE. SUBSAMPLE/WEIGHT NOTE: NHANES total cholesterol/triglyceride/LDL are drawn from the fasting sub-sample, which has its own NHANES sub-sample weight; the derivation used the measured criteria on rows where the lipid panel was present and the single general MEC exam weight WTMEC2YR for all design-weighted analyses. Because the PRIMARY analysis is unweighted complete-case, this mismatch does NOT bias the reported coefficients; it is a caveat only for the design-weighted derivation sensitivity and for lipid prevalence, which we therefore label indicative (Supplementary Table S14).
Weights: NHANES WTMEC2YR/2 (2 cycles pooled); NHIS WTFA_SA/2 (2 years pooled) for weighted descriptives; per-set weighted uses single-year weight (scale cancels for ratios).
Design: NHANES SDMVPSU/SDMVSTRA; NHIS PPSU/PSTRAT. PSU-cluster bootstrap resamples (year x stratum x pseudo-PSU).
Excluded (not in both): eyeglasses (NHANES VIQ017 is 'blind in both eyes', not eyeglasses), recent eye exam (NHIS AVISEXAM only), education/race/income (different schemes / separate NHIS personsx file).

## Construct non-equivalence (READ THIS before interpreting transportability)
Three predictors are harmonised across a **measurement-mode change** and are therefore **related but NOT strictly equivalent constructs**, not the same exposure measured twice:
- Obesity: NHANES measured BMI (objective) <-> NHIS self-reported BMI. Self-report systematically under-estimates BMI and misclassifies near the >=30 cut-point, so the NHIS "obesity" flag has lower specificity for true adiposity.
- Hypertension: NHANES measured mean SBP>=140/DBP>=90 (a physiological threshold, includes undiagnosed cases) <-> NHIS "ever told you had hypertension" (a healthcare-contact / diagnosed construct). Someone with untreated high BP is + in NHANES but - in NHIS; the two labels share a name but not a definition.
- High cholesterol: NHANES measured lipid-panel thresholds (TC>=240 / HDL<40-50 / TG>=200 / LDL>=160) <-> NHIS "ever told your cholesterol was high" (diagnosed/awareness construct), again capturing healthcare contact rather than the biochemical level.
Consequence: a coefficient estimated on the NHANES measured definition is transported onto a differently-defined NHIS variable. This mode mismatch is a *structural* source of calibration drift and attenuation of the added-value of the six clinical variables over age, independent of any real cross-population change. We therefore treat obesity/hypertension/high-cholesterol in the external set as **diagnosis-based proxies for the NHANES measurement-based exposures**, and do NOT claim the predictors are identically measured in the two surveys. Diabetes and smoking are the closest to true equivalence (both ever-diagnosed / self-reported categories on both sides).
