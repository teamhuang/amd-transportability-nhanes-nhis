#!/usr/bin/env python3
"""Phase-4 audit, ROUND-4: confirm removed overclaims/generalisations are gone, new round-4
numbers are present, abstract<=200, display items ordered, placeholders intact."""
import re, sys, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; MS=ROOT/"manuscript"/"题目4手稿.md"
txt=MS.read_text(encoding="utf-8"); low=txt.lower()

MUST_NOT=[
 "rarely report external calibration",
 "complete-case selection is mild",
 "is mild but present",
 "Transportability and calibration of a low-cost clinical risk model",  # old title
 "Discrimination is portable",
 "Measurement-mode did not explain the drift",
 "The score reflects objectively graded",
 "top-coded at 70",
 "confirmed clinically undiagnosed",
 "ELSA/SHARE/SAGE) design was originally envisioned",
]
BANNED=["principal engine","restored calibration","(bet)","load-bearing","was protective",
        "reviewer requests","v2-SR","v1-SR","v17","Running head","Significance statement","TODO","FIXME"]
REQUIRED=["0.787","0.903","0.237","2.88%","1.24","1.21",          # weighted calibration
          "0.7746","0.0932","6.06%","3.71%","1,238",              # MI + selection
          "2.67","3.09","0.836","0.971","1,000-replicate",        # survey TAYLOR/bootstrap
          "-0.0117","-0.0187","18.6%","98.1%","44.2%","93.6%",    # objective paired + concordance
          "+0.0013","0.847",                                       # RCS vs linear
          "Supplementary Table S7","Supplementary Table S8","Supplementary Table S9",
          "Supplementary Table S10","Supplementary Table S11","Supplementary Table S12","Supplementary Table S13",
          "FigSR7","FigSR8","exploratory","events per predictor parameter","35149155","32265220"]

findings=[]
for p in MUST_NOT:
    if p.lower() in low:
        ln=txt[:low.find(p.lower())].count("\n")+1
        findings.append(("BLOCKING", f"removed-phrase still present: '{p}'", f"line ~{ln}"))
for b in BANNED:
    if b.lower() in low:
        findings.append(("BLOCKING", f"banned leftover: '{b}'", "-"))
for r in REQUIRED:
    if r not in txt:
        findings.append(("WARN", f"required term/number missing: '{r}'", "-"))

m=re.search(r"## Abstract\n(.+?)\n\nKeywords", txt, re.S); aw=len(m.group(1).split()) if m else -1
if aw>200 or aw<0: findings.append(("BLOCKING", f"abstract words {aw} (>200)","-"))
title=txt.splitlines()[0].lstrip("# ").strip(); tw=len(title.split())
if tw>20: findings.append(("WARN", f"title {tw} words (>20)","-"))
supp=re.findall(r"### Supplementary Table S(\d+)", txt)
ordered = supp==sorted(supp, key=int)
if not ordered: findings.append(("BLOCKING", f"supplementary tables out of order: {supp}","-"))
report={"abstract_words":aw,"title_words":tw,"supp_order":supp,
        "author_specific_placeholders":txt.count("[AUTHOR-SPECIFIC"),
        "findings_count":len(findings),"findings":findings}
print(json.dumps(report,indent=2,ensure_ascii=False))
blocking=[f for f in findings if f[0]=="BLOCKING"]
print("\nVERDICT:", "FAIL (blocking)" if blocking else ("PASS-with-warnings" if findings else "PASS"))
sys.exit(1 if blocking else 0)
