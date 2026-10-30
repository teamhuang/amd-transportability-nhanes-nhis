"""Emit Supplementary Table S6 (unweighted-primary vs design-weighted external calibration)
from sr_weighted_calibration.json, and append to Tables_sr.md. Keeps every number code-derived.
"""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"analysis"/"output"
r=json.load(open(OUT/"sr_weighted_calibration.json"))
P=r["pooled"]; cw=r["cross_year_update"]
def pct(x): return f"{100*x:.2f}%"
lines=[]
lines.append("### Supplementary Table S6. External calibration under two estimands (frozen full model)")
lines.append("")
lines.append("| Quantity (NHIS pooled) | Unweighted (primary) | Design-weighted |")
lines.append("| --- | --- | --- |")
lines.append(f"| AUROC | {P['unweighted']['auc']:.3f} | {P['weighted']['auc']:.3f} |")
lines.append(f"| Calibration slope | {P['unweighted']['slope']:.3f} | {P['weighted']['slope']:.3f} |")
lines.append(f"| Calibration-in-the-large (logit) | {P['unweighted']['citl']:.3f} | {P['weighted']['citl']:.3f} |")
lines.append(f"| Mean predicted risk | {pct(P['unweighted']['mean_pred'])} | {pct(P['weighted']['mean_pred'])} |")
lines.append(f"| Observed prevalence | {pct(P['unweighted']['obs'])} | {pct(P['weighted']['obs'])} |")
lines.append(f"| Over-prediction ratio (pred/obs) | {P['unweighted']['overpred_ratio']:.2f} | {P['weighted']['overpred_ratio']:.2f} |")
lines.append(f"| Brier | {P['unweighted']['brier']:.4f} | {P['weighted']['brier']:.4f} |")
lines.append("")
u=cw["unweighted"]["update_2016_to_2017"]; w=cw["weighted"]["update_2016_to_2017"]
lines.append(f"Cross-year intercept update (2016→2017): unweighted mean predicted {pct(u['frozen_mean_pred'])}→{pct(u['updated_mean_pred'])} against observed {pct(u['observed'])}; design-weighted {pct(w['frozen_mean_pred'])}→{pct(w['updated_mean_pred'])} against observed {pct(w['observed'])}. Under both estimands the update moves predicted risk closer to, but below, the observed value (Brier {u['frozen_brier']:.5f}→{u['updated_brier']:.5f} unweighted; {w['frozen_brier']:.5f}→{w['updated_brier']:.5f} weighted).")
lines.append("")
wu=r["dca"]["unweighted"]; wd=r["dca"]["weighted"]
nbm_u=", ".join(f"{d['nb']:.4f}" for d in wu); nbm_w=", ".join(f"{d['nb']:.4f}" for d in wd)
nb_all_w=", ".join(f"{d['nb_treat_all']:.4f}" for d in wd)
lines.append(f"Decision-curve net benefit of the full model at thresholds 2/5/8/10%: unweighted {nbm_u}; design-weighted {nbm_w} (weighted treat-all {nb_all_w}). Net benefit remains above treat-all/treat-none at low thresholds under both estimands; its magnitude tracks the weighted prevalence and is therefore estimand-specific.")
lines.append("")
block="\n".join(lines)+"\n"
tb=OUT/"Tables_sr.md"; txt=tb.read_text(encoding="utf-8")
if "Supplementary Table S6" not in txt:
    tb.write_text(txt.rstrip()+"\n\n"+block, encoding="utf-8")
else:
    # replace from the S6 header to end
    idx=txt.find("### Supplementary Table S6")
    tb.write_text(txt[:idx].rstrip()+"\n\n"+block, encoding="utf-8")
print(block)
