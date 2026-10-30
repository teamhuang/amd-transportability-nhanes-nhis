"""ROUND-2 final: nested models, LP dump, SR tables (main+supp), and figures — single linear-age source."""
import os
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]: os.environ[v]="1"
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"; FIG=ROOT/"figures"
import numpy as np, pandas as pd, json
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, brier_score_loss
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SEED=20261028
BIN=["male","smoke_former","smoke_current","obesity","hypertension","diabetes","dyslipidemia"]
CLEAN=["male","smoke_former","smoke_current","diabetes"]
der=pd.read_csv(OUT/"nhanes_derivation.csv"); der=der[der.cc_pred==True].reset_index(drop=True)
ext=pd.read_csv(OUT/"nhis_external.csv"); ext=ext[ext.cc_pred==True].reset_index(drop=True)
def D(df,preds): return np.column_stack([np.ones(len(df)),df[['age']+preds].values.astype(float)])
def fit(preds,df): 
    X=D(df,preds); return sm.Logit(df['amd_self'].values.astype(float),X).fit(disp=0)
def metrics(y,lp):
    p=1/(1+np.exp(-lp)); y=np.asarray(y,float)
    return dict(auc=round(float(roc_auc_score(y,lp)),4),brier=round(float(brier_score_loss(y,p)),5),
        slope=round(float(sm.Logit(y,sm.add_constant(lp)).fit(disp=0).params[1]),4),
        intercept=round(float(sm.Logit(y,np.ones((len(y),1)),offset=lp).fit(disp=0).params[0]),4),
        mean_pred=round(float(p.mean()),4),obs=round(float(y.mean()),4),n=int(len(y)),events=int(y.sum()))
def ext_metrics(model,preds): 
    return metrics(ext['amd_self'].values.astype(float), D(ext,preds)@np.asarray(model.params))
nested={}
for nm,preds in [('age_only',[]),('age_sex',['male']),('mode_clean',CLEAN),('full',BIN)]:
    mo=fit(preds,der); dm=metrics(der['amd_self'].values, D(der,preds)@np.asarray(mo.params))
    nested[nm]=dict(derivation=dm,external=ext_metrics(mo,preds))
# delta AUC/Brier full vs age_only individual bootstrap (ext)
mf=fit(BIN,der); ma=fit([],der)
lf=D(ext,BIN)@np.asarray(mf.params); la=D(ext,[])@np.asarray(ma.params); ye=ext['amd_self'].values.astype(float)
rng=np.random.RandomState(SEED);dd=[]
for _ in range(1000):
    idx=rng.randint(0,len(ye),len(ye));yy=ye[idx]
    if yy.sum() in (0,len(yy)):continue
    dd.append(roc_auc_score(yy,lf[idx])-roc_auc_score(yy,la[idx]))
nested['delta_full_vs_age']=dict(dAUC=round(float(roc_auc_score(ye,lf)-roc_auc_score(ye,la)),4),
    CI=[round(float(np.percentile(dd,2.5)),4),round(float(np.percentile(dd,97.5)),4)])
json.dump(nested,open(OUT/"sr_nested.json","w"),indent=2)
# LP dump for figures
lpd=pd.DataFrame(dict(age=der['age'],lp=D(der,BIN)@np.asarray(mf.params),y=der['amd_self'],src='NHANES'))
lpe=pd.DataFrame(dict(age=ext['age'],lp=D(ext,BIN)@np.asarray(mf.params),y=ext['amd_self'],src='NHIS'))
pd.concat([lpd,lpe],ignore_index=True).to_csv(OUT/"sr_lp.csv",index=False)
json.dump(nested,open(OUT/"sr_nested.json","w"),indent=2)
print("nested:",json.dumps(nested,indent=1)[:900])

# ================= FIGURES =================
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8})
OI=dict(blue="#0072B2",orange="#E69F00",green="#009E73",verm="#D55E00",purple="#CC79A7",sky="#56B4E9",grey="#7f7f7f",black="#000000")
def save(f,n): f.savefig(FIG/(n+".png"),dpi=300,bbox_inches="tight"); f.savefig(FIG/(n+".pdf"),bbox_inches="tight"); plt.close(f)
sa=json.load(open(OUT/"sr_external_age.json"))['external']
if os.environ.get("LEGACY_FIGURES", "0") != "1":
    print("[legacy figure block skipped; figures are produced by make_figures_nature_round17.py]")
    raise SystemExit(0)

# Fig T: transportability per set
order=['NHIS_pooled','NHIS_2016','NHIS_2017','NHIS_region_Northeast','NHIS_region_Midwest','NHIS_region_South','NHIS_region_West']
names=[k.replace('NHIS_','') for k in order]; auc=[sa[k]['auc'] for k in order]; sl=[sa[k]['slope'] for k in order]; it=[sa[k]['intercept'] for k in order]
x=np.arange(len(names)); f,(a1,a2)=plt.subplots(1,2,figsize=(6.6,3.2))
a1.axhline(0.7838,ls=":",color=OI['blue'],lw=1,label="NHANES internal 0.784"); a1.bar(x,auc,0.6,color=OI['verm']); a1.set_ylim(0.5,0.9)
a1.set_xticks(x);a1.set_xticklabels(names,rotation=45,ha="right",fontsize=7);a1.set_ylabel("AUROC");a1.legend(frameon=False,fontsize=6.5)
a1.set_title("A  Discrimination transfers",fontsize=8.5)
a2.bar(x-0.2,sl,0.38,color=OI['green'],label="slope"); a2.bar(x+0.2,it,0.38,color=OI['purple'],label="intercept")
a2.axhline(1,ls="--",color=OI['black'],lw=0.7);a2.axhline(0,color=OI['grey'],lw=0.6)
a2.set_xticks(x);a2.set_xticklabels(names,rotation=45,ha="right",fontsize=7);a2.legend(frameon=False,fontsize=6.5)
a2.set_title("B  Calibration slope/intercept drift",fontsize=8.5); save(f,"FigSR1_transportability")
# Fig M: nested comparison
keys=['age_only','age_sex','mode_clean','full'];lbl=['Age','Age+sex','Mode-clean','Full']
dv=[nested[k]['derivation']['auc'] for k in keys]; ev=[nested[k]['external']['auc'] for k in keys]
f,ax=plt.subplots(figsize=(4.6,3.2));x=np.arange(4)
ax.bar(x-0.2,dv,0.38,color=OI['blue'],label="NHANES (derivation)");ax.bar(x+0.2,ev,0.38,color=OI['verm'],label="NHIS (external)")
for i in range(4): ax.text(i-0.2,dv[i]+0.004,f"{dv[i]:.3f}",ha="center",fontsize=7);ax.text(i+0.2,ev[i]+0.004,f"{ev[i]:.3f}",ha="center",fontsize=7)
ax.set_xticks(x);ax.set_xticklabels(lbl,fontsize=8);ax.set_ylim(0.6,0.85);ax.set_ylabel("AUROC");ax.legend(frameon=False,fontsize=7)
d=nested['delta_full_vs_age'];ax.set_title(f"Age alone near-optimal; dAUC(full-age)={d['dAUC']:.3f}\n(95% CI {d['CI'][0]:.3f}-{d['CI'][1]:.3f})",fontsize=8.5);save(f,"FigSR2_model_comparison")
# Fig C: calibration deciles from sr_lp
lp=pd.read_csv(OUT/"sr_lp.csv")
def dec(df):
    q=np.quantile(1/(1+np.exp(-df['lp'])),np.linspace(0,1,11));idx=np.digitize(1/(1+np.exp(-df['lp'])),q[1:-1])
    pr=[];ob=[]
    pp=1/(1+np.exp(-df['lp'].values))
    for b in range(10):
        s=idx==b
        if s.sum()>0: pr.append(pp[s].mean());ob.append(df['y'].values[s].mean())
    return pr,ob
f,ax=plt.subplots(figsize=(3.6,3.4))
for src,col,lab in [('NHANES',OI['blue'],"NHANES (derivation)"),('NHIS',OI['verm'],"NHIS (external)")]:
    pr,ob=dec(lp[lp.src==src]);ax.plot(pr,ob,"o-",color=col,label=lab)
ax.plot([0,0.12],[0,0.12],ls="--",color=OI['black'],lw=0.8,label="ideal");ax.set_xlim(0,0.12);ax.set_ylim(0,0.12)
ax.set_xlabel("Predicted probability");ax.set_ylabel("Observed AMD");ax.legend(frameon=False,fontsize=7);ax.set_title("Calibration (deciles)",fontsize=8.5);save(f,"FigSR3_calibration")
# Fig D: DCA full vs age from sr_dca
dca=json.load(open(OUT/"sr_dca.json"))
f,ax=plt.subplots(figsize=(3.8,3.2))
th=[r['threshold'] for r in dca];nb=[r['nb_model'] for r in dca];nba=[r['nb_ageonly'] for r in dca];nball=[r['nb_all'] for r in dca]
ax.plot(th,nb,"o-",color=OI['verm'],label="Full model");ax.plot(th,nba,"s--",color=OI['blue'],label="Age only")
ax.plot(th,nball,"^:",color=OI['black'],label="Treat all");ax.axhline(0,color=OI['grey'],label="Treat none")
ax.set_xlabel("Threshold probability");ax.set_ylabel("Net benefit");ax.set_title("Decision curve (NHIS)",fontsize=8.5);ax.legend(frameon=False,fontsize=6.5);save(f,"FigSR4_DCA")
# Fig O: objective
fu=json.load(open(OUT/"sr_fundus_denominators.json"))
f,(a1,a2)=plt.subplots(1,2,figsize=(6.4,3.0));lab=["any-ARM","late-AMD"]
fullv=[fu['anyARM']['auc_full'],fu['lateAMD']['auc_full']] or [0.751,0.894]
fullv=[fu['anyARM']['auc_full'],fu['lateAMD']['auc_full']];agev=[fu['anyARM']['auc_age'],fu['lateAMD']['auc_age']]
x=np.arange(2);a1.bar(x-0.2,agev,0.38,color=OI['sky'],label="Age only");a1.bar(x+0.2,fullv,0.38,color=OI['verm'],label="Full")
for i in range(2):a1.text(i-0.2,agev[i]+0.01,f"{agev[i]:.2f}",ha="center",fontsize=7);a1.text(i+0.2,fullv[i]+0.01,f"{fullv[i]:.2f}",ha="center",fontsize=7)
a1.set_xticks(x);a1.set_xticklabels(lab);a1.set_ylim(0,1);a1.axhline(0.5,color=OI['grey'],lw=0.6);a1.legend(frameon=False,fontsize=7);a1.set_title("A  Objective fundus: age >= full",fontsize=8.5)
c=fu['concordance'];a2.pie([c['fundus_pos_self_pos'],c['fundus_pos_self_neg'],c['fundus_neg_self_pos'],c['both_neg']],
    labels=["both","fundus+/report-\n(awareness gap)","report+/fundus-","neither"],colors=[OI['green'],OI['orange'],OI['verm'],OI['grey']],autopct="%1.0f",textprops={"fontsize":6.5})
a2.set_title("B  Self-report vs fundus (n=%d)"%c['n'],fontsize=8.5);save(f,"FigSR5_objective")
print("SR figures:",[x for x in sorted(os.listdir(FIG)) if x.startswith("FigSR")])
