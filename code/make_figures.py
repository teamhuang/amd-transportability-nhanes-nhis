"""Publication figures for Topic 1 (AMD clinical risk score transportability).
Arial + Okabe-Ito colorblind-safe palette; each figure saved as 300-dpi PNG + vector PDF.
All numbers are read from analysis/output (never hand-entered).
"""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[v]="1"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, json
from sklearn.metrics import roc_curve, auc, roc_auc_score

OUT=r"H:\眼科论文\题目4\analysis\output"
FIG=r"H:\眼科论文\题目4\figures"
os.makedirs(FIG,exist_ok=True)
plt.rcParams.update({
 "font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8,
 "axes.labelsize":9,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":7.5,
 "figure.dpi":150})
OI=dict(blue="#0072B2",orange="#E69F00",green="#009E73",verm="#D55E00",
        purple="#CC79A7",sky="#56B4E9",yellow="#F0E442",black="#000000",grey="#7f7f7f")

def save(fig,name):
    fig.savefig(os.path.join(FIG,name+".png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(FIG,name+".pdf"),bbox_inches="tight")
    plt.close(fig)

def load_json(n): return json.load(open(os.path.join(OUT,n)))
find=load_json("findings.json"); sens=load_json("sensitivity_objective.json")
extlp=pd.read_csv(os.path.join(OUT,"external_lp.csv"))
derp=pd.read_csv(os.path.join(OUT,"derivation_pred.csv"))
coef=pd.read_csv(os.path.join(OUT,"model_coefficients.csv"))
agedose=load_json("age_dose.json")

# ---------- Fig 1: cohort / derivation-validation flow ----------
fig,ax=plt.subplots(figsize=(5.5,3.2)); ax.axis("off")
nh=load_json("nhanes_build_summary.json"); nis=load_json("nhis_build_summary.json")
boxes=[
 (0.02,0.80,0.96,0.16,"NHANES 2005-06 + 2007-08  (derivation cohort)\nUS exam survey; objective fundus grading available",OI['blue']),
 (0.10,0.56,0.80,0.16,f"Age>=40, valid self-reported AMD:  n = {nh['n_age40_valid_amd']}\nComplete-case score model:  n = {nh['n_complete_case']}  (events = {nh['events_amd_self_cc']})",OI['sky']),
 (0.18,0.32,0.64,0.14,"Fit logistic clinical risk score (7 predictors)\nFreeze intercept + betas",OI['green']),
 (0.02,0.05,0.96,0.18,"EXTERNAL VALIDATION (fully open US survey, ever-told / self-reported)\nNHIS 2016  n=%d   +   NHIS 2017  n=%d\npooled n=%d; also per-year & per-region transportability"%(
     nis['by_year']['2016']['cc_n'],nis['by_year']['2017']['cc_n'],
     next(e['n'] for e in find['external'] if e['cohort']=='NHIS_pooled')),OI['verm'])]
for x,y,w,h,t,c in boxes:
    ax.add_patch(plt.Rectangle((x,y),w,h,fc=c,ec="black",lw=0.6,alpha=0.16))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=7.5)
for (x1,y1),(x2,y2) in [((0.5,0.80),(0.5,0.72)),((0.5,0.56),(0.5,0.46)),((0.5,0.32),(0.5,0.23))]:
    ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",lw=1.1,color="black"))
ax.set_xlim(0,1); ax.set_ylim(0,1)
save(fig,"Fig1_cohort_flow")

# ---------- Fig 2: forest plot derivation ORs (binary) + age ref ----------
lab={'male':'Male (vs female)','smoke_former':'Former smoker (vs never)','smoke_current':'Current smoker (vs never)',
     'obesity':'Obesity BMI>=30','hypertension':'Hypertension','diabetes':'Diabetes','dyslipidemia':'High cholesterol'}
sub=coef[coef['term'].isin(lab)].copy()
sub['p']=sub['p_holm']
sub=sub.sort_values('OR')
# age effect: OR at 80 vs 55 (from dose curve)
ages=agedose['ages']; ors=agedose['OR_vs_ref55']
i80=list(ages).index(80) if 80 in ages else len(ages)-1
fig,ax=plt.subplots(figsize=(5.2,3.4))
y=np.arange(len(sub))
err=np.array([np.log(sub.OR)-np.log(sub.OR_lo),np.log(sub.OR_hi)-np.log(sub.OR)])
ax.errorbar(sub.OR,y,xerr=err,fmt="o",color=OI['black'],ecolor=OI['grey'],ms=5,capsize=3,lw=1)
for k,r in enumerate(sub.itertuples()):
    star="*" if r.p_holm<0.05 else ""
    ax.text(r.OR_hi*1.25,k,f"{r.OR:.2f}{star}",va="center",fontsize=7)
ax.axvline(1,ls="--",lw=0.9,color=OI['verm'])
ax.set_yticks(y); ax.set_yticklabels([lab[t] for t in sub.term])
ax.set_xscale("log"); ax.set_xlim(0.5,4)
ax.set_xlabel("Adjusted odds ratio (95% CI), self-reported AMD")
ax.set_title("NHANES derivation: clinical predictors\n(* Holm-adjusted p<0.05)",fontsize=8.5)
ax.set_ylabel("")
# annotation for age
ax.text(0.98,0.04,f"Age effect via RCS spline: OR 80 vs 55 y = {ors[i80]:.1f}",
        transform=ax.transAxes,ha="right",fontsize=7,style="italic")
save(fig,"Fig2_forest_derivation")

# ---------- Fig 3: transportability — AUROC and calibration slope by cohort ----------
ext=find['external']
names=[e['cohort'].replace("NHIS_","") for e in ext]
aucv=[e['auc'] for e in ext]; alo=[e['auc_lo'] for e in ext]; ahi=[e['auc_hi'] for e in ext]
slp=[e['calib_slope'] for e in ext]; intr=[e['calib_intercept'] for e in ext]
x=np.arange(len(names))
fig,(a1,a2)=plt.subplots(1,2,figsize=(6.6,3.2))
a1.axhline(find['corrected_auc'],ls=":",color=OI['blue'],lw=1,label=f"Internal (NHANES) {find['corrected_auc']:.3f}")
a1.errorbar(x,aucv,yerr=[np.array(aucv)-np.array(alo),np.array(ahi)-np.array(aucv)],fmt="o",color=OI['verm'],ecolor=OI['grey'],ms=5,capsize=3,lw=1)
a1.set_xticks(x); a1.set_xticklabels(names,rotation=45,ha="right",fontsize=7)
a1.set_ylim(0.5,0.9); a1.set_ylabel("AUROC (95% CI)"); a1.legend(loc="lower right",frameon=False)
a1.set_title("A  Discrimination transfers",fontsize=8.5)
a2.axhline(1.0,ls="--",color=OI['black'],lw=0.9,label="perfect slope=1")
a2.bar(x-0.18,slp,0.34,color=OI['green'],label="Calibration slope")
a2.bar(x+0.18,intr,0.34,color=OI['purple'],label="Calib. intercept")
a2.axhline(0,color=OI['grey'],lw=0.6)
a2.set_xticks(x); a2.set_xticklabels(names,rotation=45,ha="right",fontsize=7)
a2.set_ylabel("value"); a2.legend(frameon=False,fontsize=7)
a2.set_title("B  Calibration drift (slope<1, intercept shift)",fontsize=8.5)
save(fig,"Fig3_transportability")

# ---------- Fig 4: calibration curves (deciles) derivation vs NHIS pooled ----------
def decile_cal(pred,y,nb=10):
    q=np.quantile(pred,np.linspace(0,1,nb+1)); idx=np.digitize(pred,q[1:-1])
    m=[];o=[]
    for b in range(nb):
        sel=idx==b
        if sel.sum()>0: m.append(pred[sel].mean()); o.append(y[sel].mean())
    return np.array(m),np.array(o)
fig,ax=plt.subplots(figsize=(3.6,3.4))
m,o=decile_cal(derp['pred'].values,derp['y'].values); ax.plot(m,o,"o-",color=OI['blue'],label="NHANES (derivation)")
pp=extlp[extlp.cohort=="NHIS_pooled"]; m2,o2=decile_cal(pp['pred'].values,pp['y'].values); ax.plot(m2,o2,"s-",color=OI['verm'],label="NHIS (external)")
ax.plot([0,0.12],[0,0.12],ls="--",color=OI['black'],lw=0.8,label="ideal")
ax.set_xlim(0,0.12); ax.set_ylim(0,0.12)
ax.set_xlabel("Predicted probability"); ax.set_ylabel("Observed AMD proportion")
ax.set_title("Calibration (deciles)",fontsize=8.5); ax.legend(frameon=False,fontsize=7)
save(fig,"Fig4_calibration_curve")

# ---------- Fig 5: DCA net benefit, NHIS pooled ----------
fig,ax=plt.subplots(figsize=(3.8,3.2))
for cname,col,mark in [("NHIS_pooled",OI['verm'],"-"),("NHIS_2016",OI['blue'],"--"),("NHIS_2017",OI['green'],"-.")]:
    s=extlp[extlp.cohort==cname]; y=s['y'].values; p=s['pred'].values; n=len(y); prev=y.mean()
    thr=np.linspace(0.005,0.12,120); nb=[]
    for pt in thr:
        tp=np.sum((p>=pt)&(y==1)); fp=np.sum((p>=pt)&(y==0))
        nb.append((tp-(fp*pt/(1-pt)))/n)
    ax.plot(thr,nb,col,ls=mark,lw=1.3,label=f"Model {cname.replace('NHIS_','')}")
prev=extlp[extlp.cohort=='NHIS_pooled']['y'].mean()
thr=np.linspace(0.005,0.12,120)
ax.plot(thr,prev-(1-prev)*thr/(1-thr),color=OI['black'],lw=1,label="Treat all")
ax.axhline(0,color=OI['grey'],lw=1,label="Treat none")
ax.set_xlim(0,0.12); ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit")
ax.set_title("Decision curve (self-reported AMD)",fontsize=8.5); ax.legend(frameon=False,fontsize=6.5)
save(fig,"Fig5_DCA")

# ---------- Fig 6: objective-grading sensitivity ----------
fig,(a1,a2)=plt.subplots(1,2,figsize=(6.6,3.0))
labels=["Self-report AMD\n(derivation)","Objective\nany-ARM","Objective\nlate-ARM"]
vals=[find['corrected_auc'],sens['objective_anyARM']['auc_frozen_selfreport_score_on'],sens['objective_lateARM']['auc_frozen_selfreport_score_on']]
cols=[OI['blue'],OI['orange'],OI['green']]
a1.bar(labels,vals,color=cols,width=0.6)
a1.axhline(0.5,ls="--",color=OI['grey'],lw=0.8)
for i,v in enumerate(vals): a1.text(i,v+0.01,f"{v:.3f}",ha="center",fontsize=7.5)
a1.set_ylim(0,1); a1.set_ylabel("AUROC"); a1.set_title("A  Score tracks true lesion",fontsize=8.5)
a1.tick_params(labelsize=7)
c=sens['concordance']
seg=[c['both_yes'],c['obj_only'],c['self_only'],c['neither']]
seglab=["both","objective-only\n(undiagnosed)","self-only","neither"]
a2.pie(seg,labels=seglab,colors=[OI['green'],OI['orange'],OI['verm'],OI['grey']],autopct="%1.0f",textprops={"fontsize":7})
a2.set_title("B  Self-report vs fundus AMD\n(n=%d)"%c['n'],fontsize=8.5)
save(fig,"Fig6_objective_sensitivity")

print("Figures written:",[f for f in sorted(os.listdir(FIG))])
