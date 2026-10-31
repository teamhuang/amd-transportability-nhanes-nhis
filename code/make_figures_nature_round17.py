"""Round-17: regenerate all eight submission figures in Nature style.

Specifications
- height 5.0 cm; width 5.0 cm (single panel) or 6.5 cm (two panels)
- Arial 7 pt; axes 0.6 pt; ticks 0.6 pt; no top/right spines
- Okabe-Ito colourblind-safe palette
- no titles, no in-figure notes or value labels; only bold panel letters A/B
- PDF (vector) + PNG 600 dpi, exact canvas size (bbox_inches=None)

Outputs (same file names as the manuscript):
  FigSR1_transportability, FigSR2_model_comparison, FigSR3_calibration,
  FigSR4_DCA, FigSR5_objective, FigSR6_weighted_calibration,
  FigSR7_calibration_selfreport, FigSR8_calibration_objective
"""
import json, os
for v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[v] = "1"
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "output"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

CM = 1 / 2.54
H = 5.0 * CM
W1 = 5.0 * CM
W2 = 6.5 * CM

OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", grey="#7F7F7F", black="#000000")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6,
    "axes.labelsize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5,
    "legend.fontsize": 5.0,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "lines.linewidth": 0.8,
    "lines.markersize": 2.6,
})

def canvas(width, ncols=1, wspace=0.42, left=0.16, right=0.97, bottom=0.17, top=0.88):
    fig, axes = plt.subplots(1, ncols, figsize=(width, H), squeeze=False)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
    return fig, axes[0]

def panel(fig, ax, letter):
    p = ax.get_position()
    fig.text(p.x0 - 0.06, p.y1 + 0.02, letter, fontsize=7, fontweight="bold",
             va="bottom", ha="left")

def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=600, pil_kwargs={"optimize": True})
    plt.close(fig)

# ---------------- data ----------------
sa = json.load(open(OUT / "sr_external_age.json"))["external"]
nested = json.load(open(OUT / "sr_nested.json"))
dca = json.load(open(OUT / "sr_dca.json"))
fu = json.load(open(OUT / "sr_fundus_denominators.json"))
wc = json.load(open(OUT / "sr_weighted_calibration.json"))
lp = pd.read_csv(OUT / "sr_lp.csv")

BIN = ["male", "smoke_former", "smoke_current", "obesity", "hypertension", "diabetes", "dyslipidemia"]
der = pd.read_csv(OUT / "nhanes_derivation.csv"); der = der[der.cc_pred == True].reset_index(drop=True)
ext = pd.read_csv(OUT / "nhis_external.csv"); ext = ext[ext.cc_pred == True].reset_index(drop=True)

def D(df, preds):
    return np.column_stack([np.ones(len(df)), df[["age"] + preds].values.astype(float)])

def fitlog(df, preds):
    return sm.Logit(df["amd_self"].values.astype(float), D(df, preds)).fit(disp=0)

mf = fitlog(der, BIN)
ma = fitlog(der, [])
lf = D(ext, BIN) @ np.asarray(mf.params)
la = D(ext, []) @ np.asarray(ma.params)
a16 = float(sm.Logit(ext.loc[ext.survey_year == 2016, "amd_self"].values.astype(float),
                     np.ones((int((ext.survey_year == 2016).sum()), 1)),
                     offset=lf[(ext.survey_year == 2016).values]).fit(disp=0).params[0])
lupd = lf + a16
ye = ext["amd_self"].values.astype(float)

def wilson(k, n):
    if n == 0:
        return None, None
    p = k / n; z = 1.96; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, [max(0, c - h), min(1, c + h)]

def calib_bins(lp_vec, y, nb=20):
    p = 1 / (1 + np.exp(-lp_vec))
    q = np.quantile(p, np.linspace(0, 1, nb + 1))
    idx = np.clip(np.digitize(p, q[1:-1]), 0, nb - 1)
    mp, ob, lo, hi = [], [], [], []
    for b in range(nb):
        s = idx == b; n = int(s.sum())
        if n >= 15:
            k = int(y[s].sum()); pk, pw = wilson(k, n)
            mp.append(float(p[s].mean())); ob.append(pk); lo.append(pw[0]); hi.append(pw[1])
    return np.array(mp), np.array(ob), np.array(lo), np.array(hi)

# ---------------- FigSR1 ----------------
order = ["NHIS_pooled", "NHIS_2016", "NHIS_2017", "NHIS_region_Northeast",
         "NHIS_region_Midwest", "NHIS_region_South", "NHIS_region_West"]
labels = ["Pooled", "2016", "2017", "NE", "MW", "S", "W"]
auc = [sa[k]["auc"] for k in order]
slope = [sa[k]["slope"] for k in order]
inter = [sa[k]["intercept"] for k in order]
x = np.arange(len(order))

fig, (a1, a2) = canvas(W2, 2, wspace=0.62, left=0.14, right=0.98, bottom=0.30, top=0.84)
a1.axhline(0.7838, ls=(0, (1, 1.5)), color=OI["blue"], lw=0.8, label="NHANES internal")
a1.bar(x, auc, 0.62, color=OI["verm"], edgecolor="none")
a1.set_ylim(0.5, 0.9)
a1.set_xticks(x); a1.set_xticklabels(labels, rotation=45, ha="right")
a1.set_ylabel("AUROC")
a1.legend(frameon=False, loc="upper right", handlelength=1.2, borderaxespad=0.2)

a2.bar(x - 0.19, slope, 0.36, color=OI["green"], edgecolor="none", label="Slope")
a2.bar(x + 0.19, inter, 0.36, color=OI["purple"], edgecolor="none", label="Intercept")
a2.axhline(1, ls="--", color=OI["black"], lw=0.6)
a2.axhline(0, color=OI["grey"], lw=0.5)
a2.set_ylim(-0.4, 1.0)
a2.set_xticks(x); a2.set_xticklabels(labels, rotation=45, ha="right")
a2.set_ylabel("Calibration")
a2.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
          handlelength=1.2, columnspacing=0.8, borderaxespad=0.0)
panel(fig, a1, "A"); panel(fig, a2, "B")
save(fig, "FigSR1_transportability")

# ---------------- FigSR2 ----------------
keys = ["age_only", "age_sex", "mode_clean", "full"]
lbl = ["Age", "Age+sex", "Mode-clean", "Full"]
dv = [nested[k]["derivation"]["auc"] for k in keys]
ev = [nested[k]["external"]["auc"] for k in keys]
x = np.arange(4)
fig, (ax,) = canvas(W1, 1, left=0.24, right=0.97, bottom=0.24, top=0.87)
ax.bar(x - 0.19, dv, 0.36, color=OI["blue"], edgecolor="none", label="NHANES")
ax.bar(x + 0.19, ev, 0.36, color=OI["verm"], edgecolor="none", label="NHIS")
ax.set_xticks(x); ax.set_xticklabels(lbl, rotation=20, ha="right")
ax.set_ylim(0.6, 0.85)
ax.set_ylabel("AUROC")
ax.legend(frameon=False, loc="upper right", ncol=2, handlelength=1.2,
          columnspacing=0.9, borderaxespad=0.2)
save(fig, "FigSR2_model_comparison")

# ---------------- FigSR3 ----------------
def deciles(df):
    p = 1 / (1 + np.exp(-df["lp"].values))
    q = np.quantile(p, np.linspace(0, 1, 11))
    idx = np.digitize(p, q[1:-1])
    pr, ob = [], []
    for b in range(10):
        s = idx == b
        if s.sum() > 0:
            pr.append(p[s].mean()); ob.append(df["y"].values[s].mean())
    return pr, ob

fig, (ax,) = canvas(W1, 1, left=0.24, right=0.94, bottom=0.22, top=0.88)
for src, col, lab in [("NHANES", OI["blue"], "NHANES"), ("NHIS", OI["verm"], "NHIS")]:
    pr, ob = deciles(lp[lp.src == src])
    ax.plot(pr, ob, "o-", color=col, label=lab, ms=2.8, lw=0.9)
ax.plot([0, 0.12], [0, 0.12], ls="--", color=OI["black"], lw=0.6, label="Ideal")
ax.set_xlim(0, 0.12); ax.set_ylim(0, 0.12)
ax.set_xlabel("Predicted risk"); ax.set_ylabel("Observed AMD")
ax.legend(frameon=False, loc="upper left", handlelength=1.2, borderaxespad=0.3)
save(fig, "FigSR3_calibration")

# ---------------- FigSR4 ----------------
th = [r["threshold"] for r in dca]
nb = [r["nb_model"] for r in dca]
nba = [r["nb_ageonly"] for r in dca]
nball = [r["nb_all"] for r in dca]
fig, (ax,) = canvas(W2, 1, left=0.13, right=0.97, bottom=0.17, top=0.88)
ax.plot(th, nb, "o-", color=OI["verm"], label="Full model", ms=2.4)
ax.plot(th, nba, "s--", color=OI["blue"], label="Age only", ms=2.4)
ax.plot(th, nball, "^:", color=OI["black"], label="Treat all", ms=2.4)
ax.axhline(0, color=OI["grey"], lw=0.6, label="Treat none")
ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit")
ax.legend(frameon=False, loc="lower left", handlelength=1.6, borderaxespad=0.3)
save(fig, "FigSR4_DCA")

# ---------------- FigSR5 ----------------
lab = ["Any ARM", "Late AMD"]
fullv = [fu["anyARM"]["auc_full"], fu["lateAMD"]["auc_full"]]
agev = [fu["anyARM"]["auc_age"], fu["lateAMD"]["auc_age"]]
x = np.arange(2)
fig, (a1, a2) = canvas(W2, 2, wspace=0.55, left=0.10, right=0.97, bottom=0.30, top=0.86)
a1.bar(x - 0.19, agev, 0.36, color=OI["sky"], edgecolor="none", label="Age only")
a1.bar(x + 0.19, fullv, 0.36, color=OI["verm"], edgecolor="none", label="Full model")
a1.set_xticks(x); a1.set_xticklabels(lab)
a1.set_ylim(0.5, 1.0); a1.set_ylabel("AUROC")
a1.legend(frameon=False, loc="upper right", handlelength=1.2, borderaxespad=0.2)

c = fu["concordance"]
pos_rep = c["fundus_pos_self_pos"]; pos_not = c["fundus_pos_self_neg"]
neg_rep = c["fundus_neg_self_pos"]; neg_not = c["both_neg"]
pos_tot = pos_rep + pos_not; neg_tot = neg_rep + neg_not
pct_pos = [100 * pos_rep / pos_tot, 100 * pos_not / pos_tot]
pct_neg = [100 * neg_rep / neg_tot, 100 * neg_not / neg_tot]
a2.bar([0, 1], [pct_pos[0], pct_neg[0]], 0.55, color=OI["green"], edgecolor="none", label="Reported")
a2.bar([0, 1], [pct_pos[1], pct_neg[1]], 0.55, bottom=[pct_pos[0], pct_neg[0]],
       color=OI["grey"], edgecolor="none", label="Not reported")
a2.set_xticks([0, 1]); a2.set_xticklabels(["Fundus +", "Fundus −"])
a2.set_ylim(0, 100); a2.set_ylabel("Participants (%)")
a2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1,
          handlelength=1.2, labelspacing=0.2, borderaxespad=0.0)
panel(fig, a1, "A"); panel(fig, a2, "B")
save(fig, "FigSR5_objective")

# ---------------- FigSR6 ----------------
sets = ["pooled", "NHIS_2016", "NHIS_2017"]
x = np.arange(3)
fig, (a1, a2) = canvas(W2, 2, wspace=0.55, left=0.17, right=0.97, bottom=0.18, top=0.88)
a1.plot([wc[s]["unweighted"]["obs"] for s in sets], [wc[s]["unweighted"]["mean_pred"] for s in sets],
        "o", color=OI["verm"], label="Unweighted", ms=3)
a1.plot([wc[s]["weighted"]["obs"] for s in sets], [wc[s]["weighted"]["mean_pred"] for s in sets],
        "s", color=OI["purple"], label="Weighted", ms=3)
a1.plot([0.025, 0.048], [0.025, 0.048], "--", color=OI["black"], lw=0.6)
a1.set_xlim(0.025, 0.045); a1.set_ylim(0.025, 0.048)
a1.set_xlabel("Observed risk"); a1.set_ylabel("Predicted risk")
a1.legend(frameon=False, loc="lower right", handlelength=1.0,
          labelspacing=0.2, borderaxespad=0.2)

a2.bar(x - 0.19, [wc[s]["unweighted"]["citl"] for s in sets], 0.36, color=OI["verm"], edgecolor="none", label="Unweighted")
a2.bar(x + 0.19, [wc[s]["weighted"]["citl"] for s in sets], 0.36, color=OI["purple"], edgecolor="none", label="Weighted")
a2.axhline(0, color=OI["black"], lw=0.6)
a2.set_xticks(x); a2.set_xticklabels(["Pooled", "2016", "2017"])
a2.set_ylabel("CITL")
a2.set_ylim(-0.32, 0.06)
a2.legend(frameon=False, loc="upper right", handlelength=1.2,
          labelspacing=0.2, borderaxespad=0.2)
panel(fig, a1, "A"); panel(fig, a2, "B")
save(fig, "FigSR6_weighted_calibration")

# ---------------- FigSR7 ----------------
fig, (a1, a2) = canvas(W2, 2, wspace=0.55, left=0.17, right=0.97, bottom=0.18, top=0.88)
for lp_vec, col, lab in [(lf, OI["verm"], "Full"), (la, OI["blue"], "Age only"), (lupd, OI["green"], "Updated")]:
    mp, ob, lo, hi = calib_bins(lp_vec, ye)
    a1.errorbar(mp, ob, yerr=[ob - lo, hi - ob], fmt="o-", color=col, ms=2.2, lw=0.8,
                capsize=1.0, elinewidth=0.6, label=lab)
a1.plot([0, 0.2], [0, 0.2], "--", color=OI["black"], lw=0.6, label="Ideal")
a1.set_xlim(0, 0.2); a1.set_ylim(0, 0.2)
a1.set_xlabel("Predicted risk"); a1.set_ylabel("Observed AMD")
a1.legend(frameon=False, loc="upper left", handlelength=1.4, borderaxespad=0.2)

a2.hist(1 / (1 + np.exp(-lf)), bins=40, color=OI["verm"], edgecolor="none")
a2.set_xlabel("Predicted risk"); a2.set_ylabel("Participants")
panel(fig, a1, "A"); panel(fig, a2, "B")
save(fig, "FigSR7_calibration_selfreport")

# ---------------- FigSR8 ----------------
fig, (a1, a2) = canvas(W2, 2, wspace=0.55, left=0.17, right=0.97, bottom=0.18, top=0.88)
for ax, endpoint, col in [(a1, "obj_anyARM", OI["orange"]), (a2, "obj_lateARM", OI["purple"])]:
    dd = der[der[endpoint].notna()].reset_index(drop=True)
    lp_vec = D(dd, BIN) @ np.asarray(mf.params)
    yd = dd[endpoint].values.astype(float)
    mp, ob, lo, hi = calib_bins(lp_vec, yd, nb=10)
    ax.errorbar(mp, ob, yerr=[ob - lo, hi - ob], fmt="o-", color=col, ms=2.5, lw=0.9,
                capsize=1.2, elinewidth=0.6, label="Full model")
    mx = max(np.max(mp), np.max(ob)) * 1.25
    ax.plot([0, mx], [0, mx], "--", color=OI["black"], lw=0.6, label="Ideal")
    ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_xlabel("Predicted risk"); ax.set_ylabel("Observed")
    ax.legend(frameon=False, loc="upper left", handlelength=1.4, borderaxespad=0.2)
panel(fig, a1, "A"); panel(fig, a2, "B")
save(fig, "FigSR8_calibration_objective")


# ---------------- FigSR9: cross-year recalibration ----------------
tempo = json.load(open(OUT / "sr_temporal.json"))
rows = [
    ("2017"+chr(10)+"from 2016", tempo["y2017_frozen"], tempo["y2017_recal2016"]),
    ("2016"+chr(10)+"from 2017", tempo["y2016_frozen"], tempo["y2016_recal2017"]),
]
fig, (ax,) = canvas(W2, 1, left=0.13, right=0.97, bottom=0.25, top=0.84)
for i, (lab, frozen, updated) in enumerate(rows):
    yf = frozen["mean_pred"] * 100; yu = updated["mean_pred"] * 100; yo = updated["obs"] * 100
    ax.annotate("", xy=(i, yu), xytext=(i - 0.14, yf),
                arrowprops=dict(arrowstyle="-|>", color=OI["grey"], lw=0.7, shrinkA=0, shrinkB=2))
    ax.plot(i - 0.14, yf, "s", mfc="white", mec=OI["grey"], ms=4, mew=0.8, label="Frozen" if i == 0 else None)
    ax.plot(i, yu, "o", color=OI["green"], ms=4, label="Updated" if i == 0 else None)
    ax.plot(i + 0.14, yo, "D", color=OI["blue"], ms=3.2, label="Observed" if i == 0 else None)
ax.set_xticks([0, 1]); ax.set_xticklabels([r[0] for r in rows], fontsize=5.5)
ax.set_ylim(3.2, 4.6); ax.set_ylabel("Risk (%)")
ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3,
          handlelength=1.0, columnspacing=1.2, borderaxespad=0.0)
save(fig, "FigSR9_crossyear_recalibration")

print("written:", [f"FigSR{i}" for i in range(1, 10)])