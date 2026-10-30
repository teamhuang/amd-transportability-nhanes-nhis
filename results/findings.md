# Findings & pre-delivery verdict — Topic 1

## Headline results (all emitted by code; independent recompute matched exactly)
- Derivation (NHANES 2005-08, complete-case n=5,773, 214 self-reported AMD, prevalence 3.71%, EPV≈30): apparent AUROC 0.790 -> bootstrap-corrected 0.782, Brier 0.0336. Dominant predictor = age (RCS); female vs male OR 1.60 (Holm p=0.023); high-cholesterol OR 1.31 and former-smoking OR 1.39 suggestive (nominal p 0.06/0.04) but NS after Holm.
- External transport (frozen coefficients -> NHIS 2016+2017): discrimination preserved — AUROC 0.775 (95% CI 0.762-0.787), per-year 0.775/0.774, per-region 0.761-0.786. Calibration drifted — slope 0.848, calibration-in-the-large intercept -0.233 (mean predicted 4.4% vs observed 3.6%). Drift dominated by intercept: intercept-only recalibration restored mean calibration (Brier 0.0334->0.0331) with AUROC unchanged; slope+intercept added ~nothing. Design-weighted derivation gave external AUROC 0.776 (conclusions unchanged).
- Construct validity (objective fundus grading, NHANES only): frozen self-report score AUROC 0.749 (any-ARM, n=5067) and 0.894 (late-ARM, 47 events). Fundus-positive-but-undiagnosed (316) far exceeded self-report-only (91) -> score tracks true macular disease, and quantifies under-diagnosis.
- Decision curve: model net benefit positive and above treat-all/treat-none across clinically plausible thresholds (0.0148/person at pooled prevalence).

## Scope honesty
Original cross-continental (NHANES+ELSA+SHARE+SAGE) design is NOT achievable with fully open AMD-containing data. Reported as US multi-survey transportability (NHANES -> NHIS 2016/2017), 3 individual datasets. BRFSS dropped (no current AMD item; 2008 module URL dead + non-representative). International validation flagged as future work, not simulated.

## 7-dimension self-check
1. Reproducibility / anti-double-rounding / traceability: independent recompute of AUROC (0.7748) and calibration slope (0.848) from stored LP matched findings.json; tables/manuscript generated from JSON (no hand-entry); figures regenerated after final data scripts, 0 stale. PASS.
2. Claim vs evidence: every conclusion maps to a sourced metric; NS effects labelled NS after Holm; honest nulls reported. No causal verbs. PASS.
3. Citations/numbering: 6 references verified via PubMed eutils (+2 dataset sources); unverified candidates dropped, not fabricated. Tables 1-5 and Figures 1-6 ordered by first citation. PASS.
4. Pipeline/leakage/domain static-check: analysis unit = person (no pseudo-replication); self-report outcome not contaminated by care-access predictors (eyeglasses/eye-exam excluded since absent in derivation cohort); objective grading kept as separate endpoint, never pooled; binary clinical predictors ensure transportability. PASS.
5. Statistical selection: logistic w/ RCS; Holm across predictors; ORs + CIs; bootstrap AUC CIs (1,000) + internal optimism (500) with fixed seed; Van Houwelingen slope + offset intercept; DCA; weighted sensitivity. PASS.
6. Execution hygiene / final check: BLAS/OMP threads locked to 1 before numpy; seed 20261028; deterministic reruns reproduced identical numbers; atomic writes; ran in clean re-exec. PASS.
7. Standard-readout / decision log / bet: design frozen pre-run with decision log (2013/2015 fixed-width excluded for unverifiable harmonization; NHIS 2016-17 chosen; predictor set = 7-item intersection); conclusion stated with direction + confidence + the one caveat most likely to overturn. PASS.

## Verdict
**PARTIAL PASS.** The core deliverable (derivation + internal validation + external transportability with discrimination/calibration separation + DCA + objective-grading construct validity) is complete, sourced, and reproducible. Items downgraded vs the original pitch, honestly recorded (not fabricated): (a) external validation is US-only (NHIS 2016/2017), not cross-continental; (b) BRFSS auxiliary cohort dropped; (c) education/race sensitivity and multiple-imputation not yet run; (d) author/funding/ethics/DOI placeholders pending.
