# Admission Record — Topic 1 (RUN-FREEZE, written before any model fit)

**Working title**: A no-imaging low-cost clinical risk score for self-reported age-related macular degeneration (AMD): derivation in NHANES (USA) and international-style transportability / calibration-drift across independent open population surveys.

**Honest reframing (data reality)**: The original pitch was a cross-continental validation across NHANES + ELSA + SHARE + SAGE. Under the user's explicit constraint to use ONLY fully-open, anonymous direct-download individual-level data with a harmonizable AMD self-report item, ELSA/SHARE/SAGE/HRS/CHARLS/LASI are excluded (all need registration or a DUA, and AMD self-report is unconfirmed/absent in most). BRFSS's current core questionnaire has NO AMD item (only in a 2006-2009 optional module fielded by a minority of states). The ONLY other fully-open cohort with a genuine AMD self-report is NHIS (`VIM_MDEV`, 2000-2017). Therefore the external-validation design is re-scoped to US open surveys, and the novelty is re-stated as **transportability across independent survey designs / measurement modes / eras**, NOT "across continents." This is recorded up front and must not be overclaimed in the manuscript.

**Freeze date**: 2026-10-28. After running, any change to the claim / metric / direction / tolerance / data / split / definitions is reported as a NEW analysis with the original retained; do not silently re-tune to approach a target.

## Analysis unit
Person. One record per sampled adult. No pseudo-replication (single respondent-level outcome). Survey design (weight, PSU, stratum) carried for population-representative descriptives and as a sensitivity for model/calibration; the primary prediction-model fit is unweighted (TRIPOD convention: derive individual risk; weights mainly shift intercept/variance). Design-based weighted replication is reported as a decision-log fork.

## Outcome (frozen)
**Y = self-reported AMD**: "ever told by an eye doctor / health professional that you have age-related macular degeneration (macular degeneration)."
- NHANES: `VIQ310 == 1` (cycles 2005-06, 2007-08; asked of 40+).
- NHIS: `VIM_MDEV == 1` (7/9 -> missing).
- BRFSS vision module (auxiliary): the "told you have AMD" item (exact code to be confirmed; if absent/degraded -> dropped, reported honestly).
- Do NOT pool NHANES objective fundus grading with self-report as one outcome. Objective grading is a SEPARATE sensitivity endpoint (see below).

## Predictors (frozen set; all low-cost, no imaging, no genomics, binary/continuous, transportable)
1. Age (years; restricted cubic spline, knots at ~55/65/75).
2. Sex (female reference).
3. Smoking: never (ref) / former / current (ever >=100 cigarettes for former/current split).
4. Obesity: BMI >= 30. NHANES = measured (BMXBMI); NHIS = self-reported (BMI/1000). **Measurement-mode difference is intentional and is a source of drift to quantify.**
5. Hypertension: NHANES = measured mean SBP>=140 or DBP>=90; NHIS = ever-told (`HYPEV==1`). (NHANES 2005-08 has no self-report HTN item; measurement-based.)
6. Diabetes ever-diagnosed: NHANES `DIQ010==1`; NHIS `DIBEV1==1` (borderline=3 -> not diabetic, sensitivity with borderline excluded vs as-no).
7. Dyslipidemia / high cholesterol: NHANES = measured (TC>=240 or HDL<40(M)/<50(F) or TG>=200 or LDL>=160); NHIS = ever-told (`CHLEV==1`). Measurement-mode difference flagged.
8. Regular eyeglasses / contact-lens use: NHANES `VIQ017==1`; NHIS `VIMGLASS==1`.
9. Recent eye-care visit (care-access/awareness proxy): NHANES eye-exam-in-last-12-months item (to be confirmed against VIQ dictionary, candidate VIQ051/VIQ056); NHIS `AVISEXAM` within 12 months. Included deliberately because it is the item most likely to load on care-seeking behaviour; its contribution is interrogated by the objective-grading sensitivity.

Education, race/ethnicity, income: available but EXCLUDED from the primary transportable score (cross-cohort category schemes differ and race is not a clinically actionable low-cost predictor); used only in sensitivity models, reported separately.

## Analytic sample / eligibility
- Age >= 40 (AMD clinically meaningful; NHANES VIQ/OPX asked 40+).
- Valid outcome.
- NHANES cycles 2005-06 + 2007-08 (only cycles with both AMD self-report AND objective retinal grading).
- Primary complete-case on the 7 score predictors; missingness tabled; MICE sensitivity if time permits.
- NHIS years pooled = **2016 + 2017** (clean anonymous-CSV public-use files; every item's coding verified against the actual file's value distribution before use). Each year ALSO analysed separately for era-based calibration drift, plus US-Census-region subgroups. **2013/2015 excluded (decision log, made pre-run):** FTP ships only fixed-width `.dat` with NO SAS/xpt/SPSS layout, and the AMD item name (`VIM_MDEV`->absent in 2015 layout) plus diabetes/smoking names drift year-to-year; reverse-engineering per-year fixed-width colspecs from a 384-page codebook PDF could not be verified to the integrity standard ("check not read"), so they are omitted rather than risk silent mis-mapping. 2014 also excluded (NHIS split-questionnaire redesign). Net result = 3 individual-level datasets (NHANES derivation + NHIS 2016 + NHIS 2017), satisfying the ">=2 individual datasets" requirement.

## Decision log
- **Derivation cohort = NHANES (not NHIS)**: chosen because NHANES is the only cohort with OBJECTIVE fundus grading (OPXRET) enabling the pathology-vs-care-seeking sensitivity, plus measured labs; trade-off = NHANES 2005-08 lacks eyeglasses/eye-exam-access items, so those are NOT in the score.
- **Score predictor set reduced to the 7-item intersection available in BOTH NHANES and NHIS** (age, sex, smoking, obesity, hypertension, diabetes, dyslipidemia). Eyeglasses-use (`VIMGLASS`) and recent-eye-exam (`AVISEXAM`) exist only in NHIS and were DROPPED from the score (a coefficient cannot be estimated for a predictor absent in the derivation cohort). VIQ017 is "blind in both eyes", NOT eyeglasses (label-verified) — corrected from an earlier assumption.
- **Unweighted derivation** (TRIPOD individual-risk convention); design-weighted descriptive stats and a weighted replication reported as sensitivity.
- **Education / race / income** excluded from the primary transportable score (cross-cohort category schemes differ; NHIS education lives in a separate `personsx` file) -> sensitivity only.

## Model & derivation
- Logistic regression, logit link, 7 predictors above. Derivation in NHANES. Coefficients (intercept + betas) FROZEN after derivation.
- Multiple comparisons: Holm across the 7 predictor Wald tests in the derivation model; raw p (2 sig-figs) alongside Holm-adjusted; never p=0.000.
- Effect sizes as OR with 95% CI (profile-likelihood or Wald). Age spline reported as global test + a reference OR.

## Internal validation (NHANES)
- Optimism-corrected AUROC and calibration slope via Harrell bootstrap (1,000 resamples, seed = 20261028).
- Apparent vs corrected Brier.

## External validation & transportability (frozen coefficients -> each cohort)
For NHIS-pooled, each NHIS year, and BRFSS (if used):
- Discrimination: AUROC + 95% CI (bootstrap; report method + resamples + seed).
- Calibration-in-the-large: intercept (mean predicted vs observed).
- Calibration slope: coefficient of the frozen linear predictor when refit logistic(Y ~ 1* + LP_fixed) (Van Houwelingen).
- Recalibration tiers: (a) as-is, (b) intercept-only recalibration, (c) intercept+slope recalibration -> quantify how much drift is calibration-in-the-large vs slope.
- DCA: net benefit across threshold range anchored to each cohort's AMD prevalence; compare model vs treat-all vs treat-none.
- Calibration plots (predicted vs observed by LP decile/quintile).

## Sensitivity analyses
1. **Objective-grading disambiguation (NHANES only, the intellectual core)**: evaluate the SAME frozen score against OBJECTIVE fundus AMD (OPXRET: any age-related maculopathy vs late AMD using OPDDARM/OPDDLARM/OPDDEXU/OPDDAGA) instead of self-report. If discrimination is maintained/improved against pathology (not just a care-access model), the score tracks true lesion. Report AUROC(self-report) vs AUROC(objective) and the role of the eye-care-visit predictor (re-fit without it).
2. Add education / income / race to derivation, re-check transportability.
3. Design-weighted derivation vs unweighted.
4. Exclude "borderline" diabetes; exclude the eye-care-visit predictor.

## Direction & expected effects (tolerance for PASS band)
- Derivation: age (positive, strongest), obesity/HTN/diabetes/dyslipidemia (positive, small-to-moderate OR), smoking (positive), eyeglasses & eye-care visit (positive, awareness/sed of care). AUROC in NHANES expected 0.68-0.78 (self-report).
- External validation (NHIS): AUROC expected to be preserved-ish (0.65-0.75) but calibration intercept/slope expected to DRIFT (case-mix + prevalence + measurement-mode), i.e., slope likely <1 and intercept shifting. This drift IS the finding, not a failure.
- Objective-grading sensitivity: AUROC vs fundus AMD expected LOWER than vs self-report (self-report inflated by care-seeking). Honest reporting required either way.

## Reporting
- TRIPOD (+ TRIPOD-AI not needed; state no AI/ML beyond logistic). STROBE for the observational validation.
- Effect sizes + CI only. No causal verbs ("causes/drives/leads to"); "associated with", "predicted", "transported to".
- Figures: Fig1 cohort flow/derivation-validation map; Fig2 forest of derivation ORs (Holm); Fig3 AUROC + calibration (slope/intercept) per cohort (transportability panel); Fig4 DCA; Fig5 objective-grading sensitivity.
- Every quantitative claim must be sourced to a script/JSON/log in workspace; figures regenerated by current code before delivery.

## Non-negotiable integrity rules
- Do NOT fabricate any cohort, BRFSS, or international number. A cohort is analysed only if its real individual file is downloaded and validated. If BRFSS module unusable -> drop + say so.
- Do NOT describe the design as cross-continental/international anymore; describe as multi-survey US transportability, honestly.
- Do NOT overwrite this file with a results-fit version. Findings -> `analysis/output/findings.md`.
- Do NOT pool objective grading with self-report.

## Verdict thresholds
- PASS: derivation model fits with expected dominant effects; frozen score transported to >=1 fully-open independent external cohort (NHIS) with sourced AUROC + calibration slope/intercept + DCA; era/subgroup calibration-drift quantified; objective-grading sensitivity completed in NHANES. All numbers traceable.
- PARTIAL: external validation runs but some planned cohort/year/BRFSS component unusable -> report which are done vs dropped, no fabrication.
- FAIL only if pipeline broken or numbers unsourcable.
