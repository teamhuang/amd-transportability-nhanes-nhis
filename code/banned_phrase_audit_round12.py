"""Round-12 Phase-4 banned-phrase / required-content audit for 题目4.
Scans the manuscript (which embeds the tables) plus the two generated table files.
Exits non-zero if a banned phrase is present or a required round-12 claim is missing.
"""
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = ["manuscript/题目4手稿.md", "manuscript/Supplementary_Information.md"]

BANNED = [
    "heterogeneity survived", "load-bearing", "we follow a result", "the honest conclusion",
    "reviewer requests", "reviewer-response", "the reviewer", "reviewer flagged",
    "extreme bounds", "comparable event rate", "isolates the signal", "was protective",
    "remained significant across all robustness", "actual p values are reported",
    "fully-converged single-chain", "running head", "significance statement",
    "main text excludes", "removed for sr", "s-m1", "s-m5", "s-m9", "(14 words)",
    "identity-link", "v7-sr", "v8-sr", "v9-sr", "TODO", "FIXME",
    "the codebook re-check is pending", "re-check against the live cdc codebook is pending",
    "proper mice", "proper-mice", "prespecified thresholds",
    "non-differential misclassification", "leaves the headline unchanged",
    "not driving any conclusion",
]
WHITELIST = ["version-specific doi is preferred over the concept doi"]

REQUIRED = [
    "photographic late AMD",
    "absolute predicted risks did not transfer",
    "-0.0003 to 0.0066",
    "NCHS Research Ethics Review Board",
    "QoderWork",
    "cell-by-cell",
    "fragile and not distinguishable from zero",
    "no equivalence statement",
    "Why respondents are classified on different numbers of indicators",
    "verified against the live CDC codebook",
]
ABSTRACT_MAX = 200

def words(s):
    return len([w for w in re.split(r"\s+", s) if w.strip() and not set(w) <= set("|-")])

problems = 0
combined = ""
for f in FILES:
    p = ROOT / f
    if not p.exists():
        print(f"[MISSING] {p}"); problems += 1; continue
    t = p.read_text(encoding="utf-8"); low = t.lower(); combined += "\n" + low
    print(f"\n===== {f} =====")
    for ph in BANNED:
        i = low.find(ph)
        if i == -1:
            continue
        win = low[max(0, i - 60): i + len(ph) + 60]
        legit = any(w in win for w in WHITELIST)
        if not legit:
            problems += 1
            print(f"  [BANNED] '{ph}' -> ...{win.strip()}...")

for r in REQUIRED:
    if r.lower() in combined:
        print(f"[required ok] {r}")
    else:
        print(f"[REQUIRED MISSING] {r}"); problems += 1

t = (ROOT / FILES[0]).read_text(encoding="utf-8")
m = re.search(r"## Abstract\n(.+?)\n\nKeywords", t, re.S)
if m:
    c = words(m.group(1))
    ok = c <= ABSTRACT_MAX
    problems += (0 if ok else 1)
    print(f"\n[abstract] {c} words (limit {ABSTRACT_MAX}) -> {'ok' if ok else 'TOO LONG'}")

print(f"\n{'FAIL' if problems else 'PASS'}: {problems} problem(s).")
sys.exit(1 if problems else 0)
