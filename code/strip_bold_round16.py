"""Round-16: remove bold outside titles.

Keeps the leading bold label on figure-legend lines ("- **Figure N (...).**" and
"- **Supplementary Fig. SN (...).**"); removes all other ** pairs in the manuscript and SI,
including table cells and reference text.
"""
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
files = [ROOT / "manuscript" / "题目4手稿.md", ROOT / "manuscript" / "Supplementary_Information.md"]

def strip_line(line):
    m = re.match(r"^(- \*\*(?:Figure|Supplementary Fig\.) [^\n]*?\.\*\*)(.*)$", line)
    if m:
        return m.group(1) + m.group(2).replace("**", "")
    return line.replace("**", "")

total = 0
for p in files:
    t = p.read_text(encoding="utf-8")
    n_before = t.count("**")
    new = "\n".join(strip_line(ln) for ln in t.split("\n"))
    n_after = new.count("**")
    total += n_before - n_after
    if "--apply" in sys.argv:
        p.write_text(new, encoding="utf-8")
    print(p.name, "bold markers before/after:", n_before, n_after)
print("removed markers:", total, "applied:", "--apply" in sys.argv)