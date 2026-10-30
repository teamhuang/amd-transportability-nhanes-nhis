#!/usr/bin/env python3
"""Phase-4 audit, ROUND-3: extend the round-2 audit with (i) the three reviewer-mandated
conclusion downgrades as must-NOT-appear phrases, and (ii) required design-weighted numbers.
"""
import re, sys, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; MS=ROOT/"manuscript"/"题目4手稿.md"
txt=MS.read_text(encoding="utf-8"); low=txt.lower()

# phrases that must NOT survive after the round-3 downgrades
MUST_NOT = [
    "Discrimination is portable",
    "Discrimination transports",
    "Measurement-mode did not explain the drift",
    "The score reflects objectively graded",
    "transports its discrimination across independent population surveys",
    "top-coded at 70",
    "confirmed clinically undiagnosed",
]
# round-2 banned leftovers still prohibited
BANNED = ["principal engine","restored calibration","What would change our conclusion","(bet",
          "A note on scope and honesty","load-bearing","isolates the signal","was protective",
          "reviewer requests","v2-SR","v1-SR","v17","Running head","Significance statement","TODO","FIXME"]
# canonical round-3 design-weighted numbers that MUST appear
REQUIRED_W = ["0.787","0.903","0.237","3.57%","2.88%","1.24","1.21","Supplementary Table S6","FigSR6","estimand"]
# canonical round-2 numbers that must still appear
REQUIRED_KEEP = ["0.791","0.784","0.774","0.760","0.788","0.004","522","0.856","0.894","26.8","TRIPOD","0.0596","0.0117"]

findings=[]
for p in MUST_NOT:
    if p.lower() in low:
        ln=txt[:low.find(p.lower())].count("\n")+1
        findings.append(("BLOCKING", f"overclaim phrase retained: '{p}'", f"line ~{ln}"))
for b in BANNED:
    if b.lower() in low:
        ln=txt[:low.find(b.lower())].count("\n")+1
        findings.append(("BLOCKING", f"banned leftover: '{b}'", f"line ~{ln}"))
for r in REQUIRED_W+REQUIRED_KEEP:
    if r not in txt:
        findings.append(("WARN", f"required term/number missing: '{r}'", "-"))

m=re.search(r"## Abstract\n(.+?)\n\nKeywords", txt, re.S)
abs_words=len(m.group(1).split()) if m else -1
if abs_words>200 or abs_words<0:
    findings.append(("BLOCKING", f"abstract words {abs_words} (>200 or missing)","-"))
title=txt.splitlines()[0].lstrip("# ").strip(); tw=len(title.split())
if tw>20: findings.append(("WARN", f"title {tw} words (>20)","-"))
kw=re.search(r"Keywords: (.+)", txt); kn=len([x.strip() for x in kw.group(1).split(";")]) if kw else 0
if kn>6: findings.append(("BLOCKING", f"keywords {kn} (>6)","-"))
# display items: main figures + main tables should stay modest (SR guidance ~<=8 combined is tight; flag >8 main figures)
nfig=len(set(re.findall(r"Figure (\d) ", txt)))
ntab=len(set(re.findall(r"\bTable (\d)\b", txt)))
report={"abstract_words":abs_words,"title_words":tw,"keywords":kn,
        "main_figure_refs":sorted(set(int(x) for x in re.findall(r"Figure (\d)", txt))),
        "main_table_refs":sorted(set(int(x) for x in re.findall(r"\bTable (\d)\b", txt))),
        "author_specific_placeholders":txt.count("[AUTHOR-SPECIFIC"),
        "findings_count":len(findings),"findings":findings}
print(json.dumps(report,indent=2,ensure_ascii=False))
blocking=[f for f in findings if f[0]=="BLOCKING"]
print("\nVERDICT:", "FAIL (blocking)" if blocking else ("PASS-with-warnings" if findings else "PASS"))
sys.exit(1 if blocking else 0)
