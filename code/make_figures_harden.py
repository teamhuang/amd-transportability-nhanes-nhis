"""Additional figures from hardened results (model comparison, temporal recalibration, objective vs age-only)."""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS"]: os.environ[v]="1"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd, json
OUT=r"H:\眼科论文\题目4\analysis\output"; FIG=r"H:\眼科论文\题目4\figures"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8})
OI=dict(blue="#0072B2",orange="#E69F00",green="#009E73",verm="#D55E00",purple="#CC79A7",sky="#56B4E9",grey="#7f7f7f")
h=json.load(open(os.path.join(OUT,"harden_results.json"))); m=h['models']
def save(f,n): f.savefig(os.path.join(FIG,n+".png"),dpi=300,bbox_inches="tight"); f.savefig(os.path.join(FIG,n+".pdf"),bbox_inches="tight"); plt.close(f)

# Fig7 nested model comparison (derivation vs external AUROC)
keys=["age_only","age_sex","full"]; lbl=["Age only","Age + sex","Full (7 vars)"]
dA=[m[k]['derivation']['apparent_auc'] for k in keys]; eA=[m[k]['external_pooled']['auc'] for k in keys]
x=np.arange(len(keys))
f,ax=plt.subplots(figsize=(4.4,3.2))
ax.bar(x-0.2,dA,0.38,color=OI['blue'],label="NHANES (derivation)")
ax.bar(x+0.2,eA,0.38,color=OI['verm'],label="NHIS (external)")
for i in range(len(keys)):
    ax.text(i-0.2,dA[i]+0.005,f"{dA[i]:.3f}",ha="center",fontsize=7)
    ax.text(i+0.2,eA[i]+0.005,f"{eA[i]:.3f}",ha="center",fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(lbl,fontsize=8); ax.set_ylim(0.6,0.85); ax.set_ylabel("AUROC")
ax.axhline(0.5,color=OI['grey'],lw=0.6)
d=m['delta_auc_full_minus_ageonly']
ax.set_title(f"Age alone is near-optimal\nDelta AUC (full-age) = {d['value']:.3f} ({d['ci'][0]:.3f}-{d['ci'][1]:.3f})",fontsize=8.5)
ax.legend(frameon=False,fontsize=7,loc="lower right")
save(f,"Fig7_model_comparison")

# Fig8 temporal recalibration across survey years
t=h['temporal']
f,ax=plt.subplots(figsize=(4.6,3.2))
grp=["2016 naive","2016 recal\n(from 2017)","2017 naive","2017 recal\n(from 2016)"]
vals=[t['naive_2016']['mean_pred'],t['recal2017_on_2016'] if False else t['recal2017_on2016']['mean_pred'],
      t['naive_2017']['mean_pred'],t['recal2016_on2017']['mean_pred']]
obs=[float(pd.read_csv(os.path.join(OUT,"nhis_external.csv")).query("survey_year==2016 and cc_pred==1")['amd_self'].mean()),
     None,
     float(pd.read_csv(os.path.join(OUT,"nhis_external.csv")).query("survey_year==2017 and cc_pred==1")['amd_self'].mean()),None]
xs=np.arange(4)
bars=ax.bar(xs,vals,color=[OI['orange'],OI['green'],OI['orange'],OI['green']],width=0.6)
ax.bar_label(bars,fmt="%.3f",fontsize=7)
# observed lines
o16=obs[0]; o17=obs[2]
ax.plot([0,1],[o16,o16],ls="--",color='black',lw=1)
ax.plot([2,3],[o17,o17],ls="--",color='black',lw=1)
ax.set_xticks(xs); ax.set_xticklabels(grp,fontsize=7); ax.set_ylabel("Mean predicted probability")
ax.set_title("Intercept recalibrated in the OTHER year\ndashed = observed prevalence",fontsize=8.5)
save(f,"Fig8_temporal_recalibration")

# Fig9 objective endpoint: full vs age-only + concordance
o=h['objective']['metrics']; c=h['objective']['concordance']
f,(a1,a2)=plt.subplots(1,2,figsize=(6.4,3.0))
labs=["any-ARM","late-AMD"]
fullv=[o['any-ARM']['auc_full'],o['late-AMD']['auc_full']]
agev=[o['any-ARM']['auc_ageonly'],o['late-AMD']['auc_ageonly']]
x=np.arange(2)
a1.bar(x-0.2,agev,0.38,color=OI['sky'],label="Age only"); a1.bar(x+0.2,fullv,0.38,color=OI['verm'],label="Full model")
for i in range(2):
    a1.text(i-0.2,agev[i]+0.01,f"{agev[i]:.2f}",ha="center",fontsize=7); a1.text(i+0.2,fullv[i]+0.01,f"{fullv[i]:.2f}",ha="center",fontsize=7)
a1.set_xticks(x); a1.set_xticklabels(labs); a1.set_ylim(0,1); a1.set_ylabel("AUROC (objective fundus)")
a1.axhline(0.5,color=OI['grey'],lw=0.6); a1.legend(frameon=False,fontsize=7)
a1.set_title("A  Fundus disease: age alone \u2265 full",fontsize=8.5)
a2.pie([c['fundus_pos_self_pos'],c['fundus_pos_self_neg'],c['fundus_neg_self_pos'],c['both_neg']],
       labels=["fundus+/report+","fundus+/report-\n(awareness gap)","fundus-/report+","neither"],
       colors=[OI['green'],OI['orange'],OI['verm'],OI['grey']],autopct="%1.0f",textprops={"fontsize":6.5})
a2.set_title("B  Self-report vs fundus (n=%d)"%c['n'],fontsize=8.5)
save(f,"Fig9_objective")
print("new figures:",[x for x in sorted(os.listdir(FIG)) if x.startswith(('Fig7','Fig8','Fig9'))])
