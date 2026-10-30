"""Round-15 journal-format split.

Scientific Reports limits main-text display items to 8 (figures and/or tables) and expects the
Supplementary Information as a separate file. We keep 3 main figures + 5 main tables (= 8) and move
Supplementary Tables S1-S14 and five figures (DCA, objective AUROC, weighted calibration, binned
self-report calibration, objective calibration) to Supplementary_Information.md.

The main tables are still emitted by make_tables_sr.py / make_round4_tables_figs.py; this script only
splits the generated block at "## Supplementary tables" and assembles the SI file.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "output"
MS = ROOT / "manuscript"
MANU = MS / "题目4手稿.md"
SI = MS / "Supplementary_Information.md"
MARKER = "<!-- TABLES inserted below by script; numbers emitted from analysis/output/sr_*.json -->"

tsr = (OUT / "Tables_sr.md").read_text(encoding="utf-8").strip()
tsr4 = (OUT / "Tables_sr4.md").read_text(encoding="utf-8").strip()

split = tsr.index("## Supplementary tables")
main_tables = tsr[:split].strip()
si_from_sr = tsr[split:].strip()

lines4 = tsr4.splitlines()
if lines4 and lines4[0].startswith("## Round-4"):
    lines4 = lines4[1:]
tsr4 = "\n".join(lines4).strip()

si_tables = si_from_sr + "\n\n" + tsr4

doc = MANU.read_text(encoding="utf-8")
head, block = doc.split(MARKER, 1)

# --- main-manuscript figure legends: keep 1-3, point to SI ---
fig_start = head.index("## Figure legends")
fig_block = head[fig_start:]
keep = []
for ln in fig_block.splitlines():
    if re.match(r"- \*\*Figure [123] ", ln):
        keep.append(ln)
new_fig_block = "## Figure legends\n" + "\n".join(keep) + \
    "\n- Supplementary figure legends for Supplementary Figs S1–S5 (decision curve, objective-endpoint " \
    "discrimination and concordance, weighted calibration, binned self-report calibration, and " \
    "objective-endpoint calibration) are provided in the Supplementary Information.\n"
head = head[:fig_start] + new_fig_block

si_legends = """## Supplementary figure legends
- **Supplementary Fig. S1 (FigSR4_DCA).** Decision-curve net benefit across the analysed thresholds (2%, 5%, 8%, 10%), full versus age-only versus treat-all/none.
- **Supplementary Fig. S2 (FigSR5_objective).** (A) Objective fundus AUROC, full versus age-only (any-ARM; late-AMD); (B) self-report × fundus cross-tabulation.
- **Supplementary Fig. S3 (FigSR6_weighted_calibration).** External calibration under two estimands. (A) Mean predicted versus observed risk for the frozen full model on pooled NHIS and each year, unweighted (primary) and design-weighted; (B) calibration-in-the-large intercept under both estimands, showing that the magnitude of the drift and its repair is target-population-specific.
- **Supplementary Fig. S4 (FigSR7_calibration_selfreport).** Binned calibration for the self-reported outcome on NHIS — 20 equal-count bins of predicted risk (not a LOESS/kernel smoother), each point the bin's mean predicted risk with the bin's observed proportion and Wilson 95% bars, against the 45° ideal line — for the frozen full model, the age-only model, and the intercept-updated model; (B) predicted-risk histogram (right-skewed).
- **Supplementary Fig. S5 (FigSR8_calibration_objective).** Calibration against objectively graded disease shown **separately per outcome** (any photographic ARM; late AMD), illustrating the absolute-risk mismatch when a self-report-trained score is applied to fundus-defined endpoints.
"""

si_doc = f"""# Supplementary Information

**Manuscript:** External validation and cross-survey miscalibration of a low-cost model for self-reported age-related macular degeneration across US population surveys.

**Authors:** Jianhuang Huang, Jing Wang.

## Contents
- Supplementary Tables S1–S14
- Supplementary Figs S1–S5 and their legends

{si_tables}

{si_legends}"""

new_doc = head + MARKER + "\n\n" + main_tables + "\n\n" + \
    "## Supplementary information\n\nSupplementary Tables S1–S14 and Supplementary Figs S1–S5, with their legends, " \
    "are provided in the separate Supplementary Information file.\n"

if "--apply" in sys.argv:
    MANU.write_text(new_doc, encoding="utf-8")
    SI.write_text(si_doc, encoding="utf-8")

print(json.dumps({
    "main_tables_lines": len(main_tables.splitlines()),
    "si_table_lines": len(si_tables.splitlines()),
    "main_figure_legends": len(keep),
    "si_figure_legends": 5,
    "applied": "--apply" in sys.argv,
    "main_doc_lines": len(new_doc.splitlines()),
    "si_doc_lines": len(si_doc.splitlines()),
}, indent=2))