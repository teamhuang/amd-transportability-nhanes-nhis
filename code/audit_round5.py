#!/usr/bin/env python3
"""Phase-4 audit, ROUND-5: confirm the specific inconsistencies the reviewer flagged are resolved."""
import re, sys, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; MS=ROOT/"manuscript"/"题目4手稿.md"
SI=ROOT/"manuscript"/"Supplementary_Information.md"
txt=MS.read_text(encoding="utf-8"); low=txt.lower()
si_txt=SI.read_text(encoding="utf-8") if SI.exists() else ""
combined=txt+"\n"+si_txt; combined_low=combined.lower()

MUST_NOT=[
 "<0.001",                       # stale RCS claim
 "not missing-at-random",        # MNAR overclaim
 "on no decision-relevant scale",# absolute sentence
 "non-reproducible",             # removed statistical rationale
 "AREDS clinical-progression score",  # mis-attributed
 "Transportability and calibration of a low-cost clinical risk model",  # old title
 "is mild but present",
 "Discrimination is portable",
]
REQUIRED=[
 "0.7750 vs 0.7737","+0.0013",           # RCS diff correct
 "15.0–22.7","97.6–98.4","36.8–51.8","92.8–94.2",  # Wilson CIs match S13
 "0.000003–0.00022",                       # unweighted dBrier CI lower bound shown at 6 dp (positive, excludes zero)
 "endpoint_code_audit","AREDS Report No. 18","16286620",
 "primary","key secondary","Exploratory",   # 3-tier objectives
 "offset(logit","calibration-in-the-large", # formulas
 "not survey-corrected",                    # objective bootstrap label
 # ROUND-8 additions
 "variable-type-appropriate","multinomial logit","structurally impossible",  # proper MICE
 "0 illegal","0.7744","Rubin's rules",      # proper MICE numbers
 "related but not strictly equivalent","plausible explanations for part",  # construct non-equivalence softened
 "not used to select the model",            # functional-form not selected on external
 "cross-sectional reported diagnosis","equivalent to a clinical diagnosis",  # endpoint reframing
 "method-specific, not one common design-based",  # S9 terminology
 "audit_survey_round8","1,185 clusters",    # survey code audit
 "agreement between two operational",        # S13 retitle
 "8757","1725-1730",                          # verified refs 12-13
]
findings=[]
for p in MUST_NOT:
    if p.lower() in combined_low:
        ln=combined[:combined_low.find(p.lower())].count("\n")+1
        findings.append(("BLOCKING", f"should-be-removed text present: '{p}'", f"line ~{ln}"))
for r in REQUIRED:
    if r not in combined:
        findings.append(("WARN", f"required fix/term missing: '{r}'", "-"))
m=re.search(r"## Abstract\n(.+?)\n\nKeywords", txt, re.S); aw=len(m.group(1).split()) if m else -1
if aw>200 or aw<0: findings.append(("BLOCKING", f"abstract words {aw}","-"))
title=txt.splitlines()[0].lstrip("# ").strip(); tw=len(title.split())
if tw>20: findings.append(("WARN", f"title {tw} words","-"))
supp=re.findall(r"### Supplementary Table S(\d+)", si_txt)
if supp!=sorted(supp,key=int): findings.append(("BLOCKING", f"supp order {supp}","-"))
# cross-check Wilson CIs in body against the S13 table (now in SI)
s13=re.search(r"Agreement: self-report positive among fundus-positive \| ([\d.]+)% \| ([\d.]+)-([\d.]+)%", si_txt)
body=re.search(r"sensitivity 18.6% \(95% CI ([\d.]+)–([\d.]+)\)", txt)
if s13 and body:
    if not (abs(float(body.group(1))-float(s13.group(2)))<0.05 and abs(float(body.group(2))-float(s13.group(3)))<0.05):
        findings.append(("BLOCKING", f"body vs S13 sensitivity CI mismatch: body {body.groups()} table {s13.group(2,3)}","-"))
print(json.dumps({"abstract_words":aw,"title_words":tw,"supp_order":supp,
                  "author_specific_placeholders":txt.count("[AUTHOR-SPECIFIC"),
                  "findings":findings},indent=2,ensure_ascii=False))
blocking=[f for f in findings if f[0]=="BLOCKING"]
print("\nVERDICT:", "FAIL (blocking)" if blocking else ("PASS-with-warnings" if findings else "PASS"))
sys.exit(1 if blocking else 0)
