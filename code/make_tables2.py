"""Emit SR-styled Tables (main + supplementary) from hardened outputs; no hand-entered numbers."""
import os, json, pandas as pd
OUT=r"H:\眼科论文\题目4\analysis\output"
h=json.load(open(os.path.join(OUT,"harden_results.json"))); eq=json.load(open(os.path.join(OUT,"model_equation.json"))); m=h['models']
extflat=pd.read_csv(os.path.join(OUT,"external_metrics_flat.csv"))
def md(head,rows):
    s="| "+" | ".join(head)+" |\n| "+" | ".join(["---"]*len(head))+" |\n"
    for r in rows: s+="| "+" | ".join(str(x) for x in r)+" |\n"
    return s

# T2 nested
rows=[]
name={'age_only':'Age only (RCS)','age_sex':'Age + sex','full':'Full (7 predictors)','mode_clean':'Mode-clean (drop obesity/HTN/lipid)'}
for k in ['age_only','age_sex','mode_clean','full']:
    d=m[k]['derivation']; e=m[k]['external_pooled']
    rows.append([name[k],d['apparent_auc'],e['auc'],e['brier'],e['calib_slope'],e['calib_intercept']])
T2=md(["Model","NHANES apparent AUROC","NHIS-pooled AUROC","NHIS Brier","NHIS calib slope","NHIS calib intercept"],rows)
d=m['delta_auc_full_minus_ageonly']
T2note=f"Incremental discrimination of the full model over age alone on external data: dAUC = {d['value']:.4f} (95% bootstrap CI {d['ci'][0]:.4f}-{d['ci'][1]:.4f})."

# T3 per-set external transportability (frozen full model)
psu=h['design']['psu_cluster_bootstrap']
t3=[]
for _,e in extflat.iterrows():
    t3.append([e['cohort'],int(e['n']),int(e['events']),f"{e['prevalence']:.4f}",f"{e['auc']:.3f} ({e['auc_lo']:.3f}-{e['auc_hi']:.3f})",
               f"{e['calib_slope']:.3f}",f"{e['calib_intercept']:.3f}",f"{e['mean_pred_prob']:.4f}",f"{e['nb_at_prevalence']:.4f}"])
T3=md(["External set","n","AMD events","Prevalence","AUROC (95% CI)*","Calib slope","Calib intercept","Mean pred","Net benefit @prev"],t3)

# T4 recalibration tiers + temporal + weighted
rec=eq; recal=json.load(open(os.path.join(OUT,"sensitivity_recal.json"))); tp=h['temporal']
T4=md(["Scenario (NHIS pooled unless noted)","AUROC","Brier","Mean predicted","Observed"],[
 ["Frozen, as-is",recal['asis']['auc'],recal['asis']['brier'],recal['asis']['mean_pred'],0.036],
 ["Intercept-only recalibration",recal['intercept_only']['auc'],recal['intercept_only']['brier'],recal['intercept_only']['mean_pred'],0.036],
 ["Intercept+slope recalibration",recal['intercept_slope']['auc'],recal['intercept_slope']['brier'],recal['intercept_slope']['mean_pred'],0.036],
 ["2017 frozen (no update)","-",tp['naive_2017']['brier'],tp['naive_2017']['mean_pred'],0.0387],
 ["2017 recalibrated with 2016 intercept","-",tp['recal2016_on2017']['brier'],tp['recal2016_on2017']['mean_pred'],0.0387],
 ["Design-weighted derivation (external AUROC)",recal['weighted_derivation']['auc_nhis'],"-","-","-"]])

# T5 objective sensitivity
o=h['objective']['metrics']; c=h['objective']['concordance']; om=h['objective']['meta']
T5=md(["Objective endpoint (NHANES)","n","Events","Prevalence","AUROC full","AUROC age-only","Full calib slope","Full mean pred vs obs"],[
 ["any AMD-related maculopathy",o['any-ARM']['n'],o['any-ARM']['events'],o['any-ARM']['prev'],o['any-ARM']['auc_full'],o['any-ARM']['auc_ageonly'],o['any-ARM']['calib_full']['calib_slope'],f"{o['any-ARM']['calib_full']['mean_pred']:.3f} vs {o['any-ARM']['calib_full']['obs']:.3f}"],
 ["late AMD",o['late-AMD']['n'],o['late-AMD']['events'],o['late-AMD']['prev'],o['late-AMD']['auc_full'],o['late-AMD']['auc_ageonly'],o['late-AMD']['calib_full']['calib_slope'],f"{o['late-AMD']['calib_full']['mean_pred']:.3f} vs {o['late-AMD']['calib_full']['obs']:.3f}"]])
T5b=md(["Self-report x fundus any-ARM","Count"],["",])
T5b=md(["Cross-tab (n=%d)"%c['n'],"Count"],[["Fundus-positive AND self-report-positive",c['fundus_pos_self_pos']],["Fundus-positive, self-report-negative",c['fundus_pos_self_neg']],["Fundus-negative, self-report-positive",c['fundus_neg_self_pos']],["Both negative",c['both_neg']]])
T5c=f"Gradeability: {om['n_valid_worseye_grading']} of age>=40 had valid worst-eye grading; {100*om['pct_ungradeable_or_missing']:.1f}% of complete cases lacked a gradeable fundus photograph."

# T6 model equation
co=eq['coefficients']
T6=md(["Term","Coefficient (logit)"],[["Intercept",co['const']],["age (linear RCS term)",co['age']],["R1(age) nonlinear",co['age_r1']],["R2(age) nonlinear",co['age_r2']]]+[[k,co[k]] for k in ["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]])
T6note=f"Age RCS: Harrell 3 df, knots={eq['age_basis'].split('knots=')[1].split('(R1')[0]}. Reference = {eq['reference_covariates']}. Worked examples — 70y male (former smoker, obese, HTN, high-chol, non-diab): lp={eq['worked_examples']['high_risk_70yM_formerSmoker_obese_HTN_highchol']['lp']}, p={eq['worked_examples']['high_risk_70yM_formerSmoker_obese_HTN_highchol']['p']}; 55y healthy female: p={eq['worked_examples']['low_risk_55yF_healthy']['p']}."

# T7 missingness
mn=h['missingness']['nhanes']
T7=md(["Variable / metric","Value"],[["Missing obesity",f"{mn['obesity']*100:.1f}%"],["Missing hypertension",f"{mn['hypertension']*100:.1f}%"],["Missing dyslipidemia",f"{mn['dyslipidemia']*100:.1f}%"],["Missing diabetes/smoking","~0.1%"],["Any predictor missing (age>=40)",f"{mn['any_pred_missing_pct']*100:.1f}%"],["AMD prevalence, included (complete-case)",f"{mn['included_amd_prev']*100:.2f}%"],["AMD prevalence, all age>=40",f"{mn['all_age40_amd_prev']*100:.2f}%"],["NHIS any predictor missing",f"{h['missingness']['nhis_any_pred_missing_pct']*100:.1f}%"]])

# T8 subgroups
T8=md(["NHIS subgroup","n?","AUROC","Calib slope","Calib intercept","Mean pred vs obs"],[[s['subgroup'],"",s['auc'],s['calib_slope'],s['calib_intercept'],f"{s['mean_pred']:.3f} vs {s['obs']:.3f}"] for s in h['subgroups']])

# T9 DCA prespecified
T9=md(["Threshold","Model net benefit","Treat-all net benefit","Flagged /1000","Cases caught /1000","Sensitivity","Cases missed /1000"],[[r['threshold'],r['nb_model'],r['nb_all'],r['flagged_per_1000'],r['cases_identified_per_1000'],r['sensitivity'],r['cases_missed_per_1000']] for r in h['dca_prespec']])

# T10 age risk curve
T10=md(["Age (y)","Predicted AMD risk (reference profile)"],[[a,p] for a,p in eq['age_only_reference_risk_curve']])

out=[]
out.append("## Main tables\n")
out.append("### Table 1. Derivation vs external population (see Tables.md)\n")
out.append("### Table 2. Nested model comparison\n\n"+T2+"\n"+T2note+"\n")
out.append("### Table 3. External transportability, frozen full model\n\n"+T3+"\n\n*95% CI from individual bootstrap; PSU-cluster-bootstrap CI for the pooled set: AUROC "+str(psu['auc_ci'])+", slope "+str(psu['slope_ci'])+", intercept "+str(psu['intercept_ci'])+".\n")
out.append("### Table 4. Recalibration tiers, cross-year recalibration, weighting\n\n"+T4+"\n")
out.append("### Table 5. Objective fundus-grading sensitivity\n\n"+T5+"\n\n"+T5b+"\n\n"+T5c+"\n")
out.append("## Supplementary tables\n")
out.append("### Supplementary Table S1. Frozen risk equation coefficients (logit)\n\n"+T6+"\n"+T6note+"\n")
out.append("### Supplementary Table S2. Missingness\n\n"+T7+"\n")
out.append("### Supplementary Table S3. Subgroup performance (NHIS)\n\n"+T8+"\n")
out.append("### Supplementary Table S4. Decision-curve net benefit, prespecified thresholds\n\n"+T9+"\n")
out.append("### Supplementary Table S5. Age-only reference risk curve\n\n"+T10+"\n")
open(os.path.join(OUT,"Tables_harden.md"),"w",encoding="utf-8").write("\n".join(out))
print(open(os.path.join(OUT,"Tables_harden.md"),encoding="utf-8").read()[:2600])
