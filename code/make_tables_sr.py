"""Round-2 SR tables (main + supplementary) from sr_*.json. Also writes harmonization dictionary."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"; MS=ROOT/"manuscript"
import json, pandas as pd, numpy as np
def L(n): return json.load(open(OUT/n))
nested=L("sr_nested.json"); ext=L("sr_external_age.json"); holm=L("sr_predictor_holm.json")
wprev=L("sr_weighted_prevalence.json"); fun=L("sr_fundus_denominators.json"); tempo=L("sr_temporal.json")
dca=L("sr_dca.json"); curves=L("sr_risk_curves.json"); eq=L("model_equation.json"); dv=L("sr_derivation.json")
sub=L("sr_subgroups.json"); miss=L("harden_results.json")['missingness']
def md(head,rows):
    s="| "+" | ".join(head)+" |\n| "+" | ".join(["---"]*len(head))+" |\n"
    for r in rows: s+="| "+" | ".join(str(x) for x in r)+" |\n"
    return s
psu=ext['external']['_pooled PSU-cluster CI']; pooled=ext['external']['NHIS_pooled']

# --- full-precision mean predicted risk on the analytic sample ---
# The stored 5-dp means are not precise enough for 2-dp percentage display: 2016 is
# 0.043354 (4.34%) but is stored as 0.04335, which floats to 4.334999 and formats as 4.33.
_co=eq['coefficients_full_precision']; _ks=['age','male','smoke_former','smoke_current','obesity','hypertension','diabetes','dyslipidemia']
_x=pd.read_csv(OUT/'nhis_external.csv'); _x=_x[_x['cc_pred']==True].copy()
_b=np.array([_co['const']]+[_co[k] for k in _ks])
_x['_p']=1/(1+np.exp(-(_b[0]+_x[_ks].values.astype(float)@_b[1:])))
_reg={1:'Northeast',2:'Midwest',3:'South',4:'West'}
def full_mean(key):
    if key=='NHIS_pooled': g=_x
    elif key.startswith('NHIS_region_'):
        code={v:k for k,v in _reg.items()}[key.replace('NHIS_region_','')]; g=_x[_x['region']==code]
    else: g=_x[_x['survey_year']==int(key.split('_')[1])]
    return float(g['_p'].mean())
# integrity check: full-precision means must agree with the stored 5-dp values
for _k in ['NHIS_pooled','NHIS_2016','NHIS_2017']+[f'NHIS_region_{v}' for v in _reg.values()]:
    _stored=ext['external'][_k]['mean_pred']; _calc=full_mean(_k)
    assert abs(_stored-_calc)<5e-5, (_k,_stored,_calc)

# T1 characteristics weighted vs unweighted
rows=[]
lab=[("Age, mean (y)","age_unw","age_w"),("Female (%)","female_unw","female_w"),("Current smoker (%)","smoke_current_unw","smoke_current_w"),
 ("Obesity BMI>=30 (%)","obesity_unw","obesity_w"),("Hypertension (%)","hypertension_unw","hypertension_w"),
 ("Diabetes (%)","diabetes_unw","diabetes_w"),("High cholesterol (%)","dyslipidemia_unw","dyslipidemia_w"),("Self-reported AMD (%)","amd_self_unw","amd_self_w")]
for name,uk,wk in lab:
    rows.append([name, wprev['NHANES'].get(uk), wprev['NHANES'].get(wk,"—"), wprev['NHIS'].get(uk), wprev['NHIS'].get(wk,"—")])
rows.insert(0,["N (>=40, complete-case)",dv['derivation_n'],"",pooled['n'],""])
T1=md(["Characteristic","NHANES unweighted","NHANES weighted","NHIS unweighted","NHIS weighted"],rows)

# T2 nested
T2=md(["Model","Derivation AUROC","NHIS-pooled AUROC","NHIS slope","NHIS intercept","NHIS Brier"],
  [[n,dv_['derivation']['auc'],dv_['external']['auc'],dv_['external']['slope'],dv_['external']['intercept'],dv_['external']['brier']]
   for n,dv_ in [("Age only",nested['age_only']),("Age+sex",nested['age_sex']),("Mode-clean",nested['mode_clean']),("Full (7)",nested['full'])]])
inc=L("sr_increment.json")
_u=inc["increment"]["unweighted"]; _w=inc["increment"]["weighted"]
T2n=(f"Incremental performance of the full model over age only (external; single frozen 1,000-replicate paired bootstrap, seed 20261028). "
     f"Unweighted analytic-sample estimand: dAUROC {_u['dAUC']:.4f} (95% CI {_u['dAUC_CI'][0]:.4f}-{_u['dAUC_CI'][1]:.4f}); "
     f"dBrier {_u['dBrier']:.6f} (95% CI {_u['dBrier_CI'][0]:.6f} to {_u['dBrier_CI'][1]:.6f}); "
     f"dNetBenefit@5% {_u['dNB5']:.5f} (95% CI {_u['dNB5_CI'][0]:.5f} to {_u['dNB5_CI'][1]:.5f}). "
     f"Design-weighted estimand: dAUROC {_w['dAUC']:.4f} (95% CI {_w['dAUC_CI'][0]:.4f}-{_w['dAUC_CI'][1]:.4f}); "
     f"dBrier {_w['dBrier']:.6f} (95% CI {_w['dBrier_CI'][0]:.6f} to {_w['dBrier_CI'][1]:.6f}). "
     f"All dAUROC/dBrier/dNB values in this manuscript come from this one run; the PSU-cluster interval for the AUROC point estimate (not the difference) is given in Table 3.")

# T3 external per-set
order=['NHIS_pooled','NHIS_2016','NHIS_2017','NHIS_region_Northeast','NHIS_region_Midwest','NHIS_region_South','NHIS_region_West']
T3=md(["External set (frozen full model)","n","Events","Prevalence (%)","AUROC (indiv. boot 95% CI)","Calib slope","Calib intercept","Brier","Mean pred vs obs (%)"],
 [[k.replace("NHIS_",""),ext['external'][k]['n'],ext['external'][k]['events'],f"{ext['external'][k]['events']/ext['external'][k]['n']*100:.2f}",
   f"{ext['external'][k]['auc']:.3f} ({ext['external'][k]['auc_ci_individual'][0]:.3f}-{ext['external'][k]['auc_ci_individual'][1]:.3f})",
   ext['external'][k]['slope'],ext['external'][k]['intercept'],ext['external'][k]['brier'],f"{full_mean(k)*100:.2f} vs {ext['external'][k]['events']/ext['external'][k]['n']*100:.2f}"] for k in order])
T3n=f"Pooled **PSU-cluster-bootstrap** 95% intervals: AUROC {psu['auc']}, slope {psu['slope']}, intercept {psu['intercept']} (composite-cluster bootstrap intervals — an approximation to design-based inference, not a Taylor/`svyglm` variance; the AUROC CI shown in the table is the individual-level bootstrap)."

# T4 recalibration/temporal/weighted (frozen rows reuse the Table-3 transport output)
T4=md(["Scenario","Brier","Mean predicted (%)","Observed (%)"],[
 ["NHIS 2017 frozen",tempo['y2017_frozen']['brier'],f"{full_mean('NHIS_2017')*100:.2f}",f"{ext['external']['NHIS_2017']['events']/ext['external']['NHIS_2017']['n']*100:.2f}"],
 ["NHIS 2017, intercept updated from 2016",tempo['y2017_recal2016']['brier'],f"{tempo['y2017_recal2016']['mean_pred']*100:.2f}",f"{tempo['y2017_recal2016']['obs']*100:.2f}"],
 ["NHIS 2016 frozen",tempo['y2016_frozen']['brier'],f"{full_mean('NHIS_2016')*100:.2f}",f"{ext['external']['NHIS_2016']['events']/ext['external']['NHIS_2016']['n']*100:.2f}"],
 ["NHIS 2016, intercept updated from 2017",tempo['y2016_recal2017']['brier'],f"{tempo['y2016_recal2017']['mean_pred']*100:.2f}",f"{tempo['y2016_recal2017']['obs']*100:.2f}"]])
T4pre="Frozen rows reuse the same transport-path output as Table 3, so the two tables agree exactly; intercept-updated rows use the temporal-update path. All risks are on the percentage scale used in the text."
T4n=T4pre+f" Intercept corrections: 2016 = {tempo['intercept_2016']}, 2017 = {tempo['intercept_2017']}; frozen-model AUROC 2016={tempo['auroc_frozen_2016']}, 2017={tempo['auroc_frozen_2017']} (unchanged by intercept shift). Design-weighted derivation external AUROC = {L('sensitivity_recal.json')['weighted_derivation']['auc_nhis']}. Note: the two cross-year directions behave differently — updating 2017 from 2016 overshoots into slight UNDER-prediction (updated {tempo['y2017_recal2016']['mean_pred']*100:.2f}% vs observed {tempo['y2017_recal2016']['obs']*100:.2f}%), while updating 2016 from 2017 still OVER-predicts (updated {tempo['y2016_recal2017']['mean_pred']*100:.2f}% vs observed {tempo['y2016_recal2017']['obs']*100:.2f}%); an intercept-only update zeroes calibration-in-the-large **within the sample used to estimate it** but never corrects the slope, so applied to an independent year it does not reliably restore calibration."

# T5 objective
T5=md(["Objective endpoint (NHANES)","n","Events","Prevalence (%)","AUROC full","AUROC age-only","Full calib slope","Full mean pred vs obs (%)"],
 [["any-ARM",fun['anyARM']['n'],fun['anyARM']['events'],f"{fun['anyARM']['events']/fun['anyARM']['n']*100:.2f}",fun['anyARM']['auc_full'],fun['anyARM']['auc_age'],fun['anyARM']['calib']['slope'],f"{fun['anyARM']['calib']['mean_pred']*100:.2f} vs {fun['anyARM']['calib']['obs']*100:.2f}"],
  ["late-AMD",fun['lateAMD']['n'],fun['lateAMD']['events'],f"{fun['lateAMD']['events']/fun['lateAMD']['n']*100:.2f}",f"{fun['lateAMD']['auc_full']} (CI {fun['lateAMD']['auc_full_CI95']})",fun['lateAMD']['auc_age'],fun['lateAMD']['calib']['slope'],f"{fun['lateAMD']['calib']['mean_pred']*100:.2f} vs {fun['lateAMD']['calib']['obs']*100:.2f}"]])
T5b=md(["Cross-tab (n=%d)"%fun['concordance']['n'],"Count"],[["Fundus+ & self-report+",fun['concordance']['fundus_pos_self_pos']],["Fundus+ , self-report- (not reporting a diagnosis)",fun['concordance']['fundus_pos_self_neg']],["Fundus- , self-report+",fun['concordance']['fundus_neg_self_pos']],["Both -",fun['concordance']['both_neg']]])
T5c=f"Denominators: among {dv['derivation_n']} analytic complete cases, {fun['n_valid_anyARM_within_cc']} had valid worst-eye any-ARM grading ({fun['pct_no_anyARM_grading_within_cc']}% ungradeable) and {fun['n_valid_lateARM_within_cc']} valid late-ARM grading; the two endpoints differ slightly because photograph quality/eye availability differ by grading variable. Late-AMD has only {fun['lateAMD']['events']} events (AUROC CI {fun['lateAMD']['auc_full_CI95']})."

# S1 equation
import math
co=eq['coefficients_4dp']
s1rows=[["Intercept",co['const'],"—"]]+[[t,co[t],round(math.exp(co[t]),3)] for t in ['age','male','smoke_former','smoke_current','obesity','hypertension','diabetes','dyslipidemia']]
S1=md(["Term","Coefficient (logit)","OR"],s1rows)
S1n=f"Age entered as a LINEAR term (OR {dv['per_year_OR']} per year, {dv['per_decade_OR']} per decade); a 3-df restricted cubic spline was inspected only post hoc (raised external AUROC by +0.0013, 0.7750 vs 0.7737) and was **not** used to select the published linear term. Reference profile: {eq['reference']}. Worked examples: a 70-year-old man who is a former smoker, obese, hypertensive and hypercholesterolaemic, without diabetes, has predicted risk {eq['worked_examples']['70yM(former,obese,HTN,highchol,nonDiab)']}; a healthy 55-year-old woman has {eq['worked_examples']['55F healthy']}. Reproduction: with 4-dp coefficients the maximum absolute deviation in linear predictor across the cohort is {eq['repro_check']['max_abs_delta_LP_4dp_vs_full']} (probability {eq['repro_check']['max_abs_delta_prob']}), so `risk_calculator.py` reproduces the model. **Locations of attached artefacts:** the full-precision coefficients for both the linear and spline specifications and the spline knots are stored in `model_equation.json` (this table shows 4-dp values for readability); the executable calculator is `risk_calculator.py`; the per-predictor joint Wald tests with Holm-adjusted p-values are in the table below (Part B of S1)."
# S1 Part B: per-predictor joint Wald test + Holm (from sr_predictor_holm.json)
def _p(x): return ("0.000" if x<=0 else (f"{x:.1e}" if x<0.001 else f"{x:.3f}"))
wald_order=[("age","Age (1 df)"),("male","Male vs female (1 df)"),("smoke","Smoking, joint former+current (2 df)"),("obesity","Obesity BMI>=30 (1 df)"),("hypertension","Hypertension (1 df)"),("diabetes","Diabetes (1 df)"),("dyslipidemia","High cholesterol (1 df)")]
S1w=md(["Predictor (joint Wald)","chi-square","df","p","Holm-adjusted p"],[[lbl,holm[k]['chi2'],holm[k]['df'],_p(holm[k]['p']),_p(holm[k]['p_holm'])] for k,lbl in wald_order])
S1w="**S1 Part B — joint Wald tests and Holm-adjusted p (Holm across the seven predictors):**\n\n"+S1w

# S2 missingness
S2=md(["Variable / metric","NHANES"],[["Missing obesity",f"{miss['nhanes']['obesity']*100:.1f}%"],["Missing hypertension",f"{miss['nhanes']['hypertension']*100:.1f}%"],["Missing dyslipidemia",f"{miss['nhanes']['dyslipidemia']*100:.1f}%"],["Missing diabetes or smoking","~0.1%"],["Any predictor missing (age>=40)",f"{miss['nhanes']['any_pred_missing_pct']*100:.1f}%"],["AMD prevalence: included vs all age>=40",f"{miss['nhanes']['included_amd_prev']*100:.2f}% vs {miss['nhanes']['all_age40_amd_prev']*100:.2f}%"],["NHIS any predictor missing",f"{miss['nhis_any_pred_missing_pct']*100:.1f}%"]])

# S3 subgroups (with n)
S3=md(["NHIS subgroup","n","Events","AUROC (95% CI)","slope","intercept","Mean pred vs obs (%)"],
 [[s['subgroup'],s['n'],s['events'],f"{s['auc']:.3f} ({s['auc_ci_individual'][0]:.3f}-{s['auc_ci_individual'][1]:.3f})",s['slope'],s['intercept'],f"{s['mean_pred']*100:.2f} vs {s['prev']*100:.2f}"] for s in sub])

# S4 DCA
S4=md(["Threshold","Flagged/1000","Cases caught/1000","Sensitivity","Cases missed/1000","NB full","NB age-only","NB all","NB none"],
 [[r['threshold'],r['flagged_per_1000'],r['cases_caught_per_1000'],r['sensitivity'],r['cases_missed_per_1000'],r['nb_model'],r['nb_ageonly'],r['nb_all'],0] for r in dca])

# S5 risk curves
S5=md(["Age","Full-model risk (reference profile)","Age-only model risk"],[[c[0],c[1],curves['age_only_model'][i][1]] for i,c in enumerate(curves['full_model_reference_profile'])])

out="## Tables\n\n### Table 1. Derivation vs external population (unweighted and design-weighted)\n\n"+T1+"\n\n### Table 2. Nested models (frozen, linear-age canonical)\n\n"+T2+"\n\n"+T2n+"\n\n### Table 3. External transportability, frozen full model\n\n"+T3+"\n\n"+T3n+"\n\n### Table 4. Recalibration / cross-year update / weighting\n\n"+T4+"\n\n"+T4n+"\n\n### Table 5. Objective fundus-grading sensitivity\n\n"+T5+"\n\n"+T5b+"\n\n"+T5c+"\n\n## Supplementary tables\n\n### Supplementary Table S1. Frozen risk equation (logit)\n\n"+S1+"\n\n"+S1n+"\n\n"+S1w+"\n\n### Supplementary Table S2. Missingness\n\n"+S2+"\n\n### Supplementary Table S3. Subgroup performance (NHIS)\n\n"+S3+"\n\n### Supplementary Table S4. Decision-curve at the analysed thresholds (2%, 5%, 8%, 10%)\n\n"+S4+"\n\n### Supplementary Table S5. Predicted risk by age\n\n"+S5+"\n"
open(OUT/"Tables_sr.md","w",encoding="utf-8").write(out)

# harmonization dictionary
harm="""# Variable harmonization dictionary (NHANES derivation <-> NHIS external validation)
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
"""
open(MS/"harmonization_dictionary.md","w",encoding="utf-8").write(harm)
print("wrote Tables_sr.md and harmonization_dictionary.md")
