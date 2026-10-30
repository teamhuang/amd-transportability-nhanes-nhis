"""Generate Markdown + CSV result tables for Topic 1 from the JSON/CSV outputs (no hand-entered numbers)."""
import os, json, numpy as np, pandas as pd
OUT=r"H:\眼科论文\题目4\analysis\output"
def J(n): return json.load(open(os.path.join(OUT,n)))
find=J("findings.json"); internal=J("internal_validation.json"); sens=J("sensitivity_objective.json")
recal=J("sensitivity_recal.json"); agedose=J("age_dose.json")
coef=pd.read_csv(os.path.join(OUT,"model_coefficients.csv"))
extflat=pd.read_csv(os.path.join(OUT,"external_metrics_flat.csv"))
der=pd.read_csv(os.path.join(OUT,"nhanes_derivation.csv")); der=der[der.cc_pred==True]
ext=pd.read_csv(os.path.join(OUT,"nhis_external.csv")); ext=ext[ext.cc_pred==True]

def md_table(header,rows):
    s="| "+" | ".join(header)+" |\n"
    s+="| "+" | ".join(["---"]*len(header))+" |\n"
    for r in rows: s+="| "+" | ".join(str(x) for x in r)+" |\n"
    return s

# Table 1 characteristics
lab=[("n (>=40, complete-case)","N"),("Age, mean (SD)","age"),("Female (%)","female"),
     ("Former smoker (%)","smoke_former"),("Current smoker (%)","smoke_current"),
     ("Obesity BMI>=30 (%)","obesity"),("Hypertension (%)","hypertension"),
     ("Diabetes (%)","diabetes"),("High cholesterol (%)","dyslipidemia"),("Self-reported AMD (%)","amd_self")]
def col_stats(df,key):
    if key=="N": return f"{len(df)}"
    if key=="female": return f"{100*(1-df['male'].mean()):.1f}"
    if key=="age": return f"{df['age'].mean():.1f} ({df['age'].std():.1f})"
    if key=="amd_self": return f"{100*df['amd_self'].mean():.2f}"
    return f"{100*df[key].mean():.1f}"
rows=[[t, col_stats(der,k), col_stats(ext,k)] for t,k in lab]
t1=md_table(["Characteristic","NHANES (derivation)","NHIS 2016-2017 (external)"],rows)

# Table 2 derivation binary predictors
binrows=[]
name={'male':'Male (ref: female)','smoke_former':'Former smoker (ref: never)','smoke_current':'Current smoker (ref: never)',
      'obesity':'Obesity BMI>=30','hypertension':'Hypertension (measured >=140/90)','diabetes':'Diabetes (ever dx)',
      'dyslipidemia':'High cholesterol / dyslipidemia'}
for p in ['male','smoke_former','smoke_current','obesity','hypertension','diabetes','dyslipidemia']:
    r=coef[coef.term==p].iloc[0]
    ph = "-" if pd.isna(r.p_holm) else (f"{r.p_holm:.3g}" if r.p_holm>=0.001 else "<0.001")
    pr = f"{r.p_raw:.3g}" if r.p_raw>=0.001 else "<0.001"
    binrows.append([name[p],f"{r.OR:.2f} ({r.OR_lo:.2f}-{r.OR_hi:.2f})",pr,ph])
# age effect row
i55=agedose['ages'].index(55); i80=agedose['ages'].index(80)
binrows.append(["Age (RCS spline); OR 80 vs 55 y",f"{agedose['OR_vs_ref55'][i80]:.1f} (vs 55)","-","dominant"])
t2=md_table(["Predictor","Adjusted OR (95% CI)","p (raw)","p (Holm)"],binrows)

# Table 3 external validation metrics
t3rows=[]
for _,e in extflat.iterrows():
    t3rows.append([e['cohort'],int(e['n']),int(e['events']),f"{e['prevalence']:.4f}",
        f"{e['auc']:.3f} ({e['auc_lo']:.3f}-{e['auc_hi']:.3f})",
        f"{e['calib_slope']:.3f}",f"{e['calib_intercept']:.3f}",f"{e['brier']:.4f}",
        f"{e['mean_pred_prob']:.4f}",f"{e['nb_at_prevalence']:.4f}"])
t3=md_table(["External set","n","AMD events","Prevalence","AUROC (95% CI)","Calib slope","Calib intercept","Brier","Mean pred prob","Net benefit @prev"],t3rows)

# Table 4 recalibration tiers + weighted
t4=md_table(["Scenario","AUROC","Brier","Mean predicted prob"],[
 ["Frozen model, as-is",recal['asis']['auc'],recal['asis']['brier'],recal['asis']['mean_pred']],
 ["Intercept-only recalibration",recal['intercept_only']['auc'],recal['intercept_only']['brier'],recal['intercept_only']['mean_pred']],
 ["Intercept+slope recalibration",recal['intercept_slope']['auc'],recal['intercept_slope']['brier'],recal['intercept_slope']['mean_pred']],
 ["Design-weighted derivation (sensitivity)",recal['weighted_derivation']['auc_nhis'],"-","-"]])

# Table 5 objective-grading sensitivity
o=sens['objective_anyARM']; l=sens['objective_lateARM']; c=sens['concordance']
t5=md_table(["Endpoint (NHANES only)","n","Events","Prevalence","AUROC of frozen self-report score (95% CI)"],[
 ["Self-reported AMD (derivation)",find['derivation_n'],find['events'],f"{find['prev']:.4f}",f"{find['corrected_auc']:.3f} (bootstrap-corrected)"],
 ["Objective any-ARM (fundus)",o['n'],o['events'],f"{o['prev']:.4f}",f"{o['auc_frozen_selfreport_score_on']:.3f} ({o['ci'][0]:.3f}-{o['ci'][1]:.3f})"],
 ["Objective late-ARM (fundus)",l['n'],l['events'],f"{l['prev']:.5f}",f"{l['auc_frozen_selfreport_score_on']:.3f} ({l['ci'][0]:.3f}-{l['ci'][1]:.3f})"]])
t5b=md_table(["Concordance self-report vs fundus any-ARM","Count"],[["Both positive",c['both_yes']],["Fundus-positive only (undiagnosed)",c['obj_only']],["Self-report-positive only",c['self_only']],["Neither",c['neither']],["Total (both measured)",c['n']]])

with open(os.path.join(OUT,"Tables.md"),"w") as f:
    f.write("### Table 1. Study population characteristics\n\n"+t1+"\n")
    f.write("### Table 2. NHANES derivation model (risk equation) — adjusted effects\n\n"+t2+"\n")
    f.write("### Table 3. External validation (transportability) performance by set\n\n"+t3+"\n")
    f.write("### Table 4. Recalibration tiers and weighted sensitivity (NHIS pooled)\n\n"+t4+"\n")
    f.write("### Table 5. Objective fundus-grading sensitivity (NHANES)\n\n"+t5+"\n\n"+t5b+"\n")
print(open(os.path.join(OUT,"Tables.md")).read())
