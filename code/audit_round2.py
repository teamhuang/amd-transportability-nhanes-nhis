#!/usr/bin/env python3
"""Phase-4 compliance + banned-phrase audit for round-2 SR manuscript.
Deterministic: reads the manuscript text from disk, flags revision leftovers,
confirms canonical required numbers appear, checks abstract word count and title.
"""
import re, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # .../题目4
MS = ROOT / "manuscript" / "题目4手稿.md"
txt = MS.read_text(encoding="utf-8")

BANNED = [
    "top-coded at 70", "principal engine", "restored calibration",
    "confirmed clinically undiagnosed", "clinically undiagnosed patients",
    "What would change our conclusion", "bet)", "(bet", "A note on scope and honesty",
    "decision-eligible", "load-bearing", "isolates the signal", "was protective",
    "The honest conclusion", "reviewer requests", "v2-SR", "v1-SR", "v17",
    "heterogeneity survived", "extreme bounds", "comparable event rate",
    "Running head", "Significance statement", "EJE", "Main text excludes",
    "we found that hypertension was protective", "10-0", "TODO", "FIXME",
]
# word-count "(NN words)" style working notes
BANNED_REGEX = [r"\(\d+\s+words?\)", r"round\s*\d+\s*audit", r"PLACEHOLDER"]

# canonical numbers that MUST appear (from sr_*.json / Tables_sr.md)
REQUIRED = [
    "0.791", "0.784", "0.774", "0.760", "0.788",      # AUROC + PSU CI
    "0.87", "0.213", "4.4", "3.6",                    # slope, intercept, mean pred/obs
    "0.004", "0.001", "0.007",                        # incremental AUROC CI
    "158.1", "2.57", "0.014",                         # joint Wald, per-decade OR, Holm male
    "522", "0.856",                                   # DCA flagged, sensitivity
    "0.894", "47", "0.849", "0.933",                  # late AMD AUROC, events, CI
    "5,067", "5,100", "316", "91",                    # fundus denominators + cross-tab
    "26.8",                                           # EPV
    "TRIPOD", "PSU", "0.0596", "0.0117",              # reporting, design-based, worked examples
]

findings = []
low = txt.lower()
for b in BANNED:
    if b.lower() in low:
        # locate
        idx = low.find(b.lower())
        line = txt[:idx].count("\n") + 1
        findings.append(("BLOCKING", f"banned phrase '{b}'", f"line ~{line}"))
for rg in BANNED_REGEX:
    for m in re.finditer(rg, txt, re.IGNORECASE):
        line = txt[:m.start()].count("\n") + 1
        findings.append(("BLOCKING", f"banned pattern '{rg}' -> '{m.group()}'", f"line ~{line}"))

missing = [r for r in REQUIRED if r not in txt]
for r in REQUIRED:
    if r not in txt:
        findings.append(("WARN", f"required number/term '{r}' not found", "-"))

# abstract word count
m = re.search(r"## Abstract\n(.+?)\n\nKeywords", txt, re.S)
abs_words = len(m.group(1).split()) if m else -1
if abs_words > 200 or abs_words < 0:
    findings.append(("BLOCKING", f"abstract word count {abs_words} (>200 or missing)", "-"))

# title word count (<=20)
title = txt.splitlines()[0].lstrip("# ").strip()
tw = len(title.split())
if tw > 20:
    findings.append(("WARN", f"title {tw} words (>20)", "-"))

# keywords count (<=6)
km = re.search(r"Keywords: (.+)", txt)
kw = [k.strip() for k in km.group(1).split(";")] if km else []
if len(kw) > 6:
    findings.append(("BLOCKING", f"keywords {len(kw)} (>6)", "-"))

# leftover AUTHOR-SPECIFIC should still exist (that's expected, report as OBSERVATION count)
n_placeholder = txt.count("[AUTHOR-SPECIFIC")

# table/figure numbering order in body
figrefs = [int(x) for x in re.findall(r"Figure (\d)", txt)]
tabrefs = [int(x) for x in re.findall(r"\bTable (\d)\b", txt)]

report = {
    "abstract_words": abs_words,
    "title_words": tw,
    "keywords": len(kw),
    "author_specific_placeholders": n_placeholder,
    "figure_first_appearance_order": figrefs[:8],
    "table_reference_order": sorted(set(tabrefs)),
    "blocking_or_warn_count": len(findings),
    "findings": findings,
}
print(json.dumps(report, indent=2, ensure_ascii=False))
blocking = [f for f in findings if f[0] == "BLOCKING"]
print("\nVERDICT:", "FAIL (blocking present)" if blocking else ("PASS-with-warnings" if findings else "PASS"))
sys.exit(1 if blocking else 0)
