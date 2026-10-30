"""Round-16 renumbering by order of first appearance.

- References 1-16 -> order of first in-text citation (refs 7/8 now cited in Methods).
- Supplementary Tables S1-S14 -> order of first citation (S2/S5 newly cited).
- Main figures 1-3: swap 2 and 3 so first citation order is 1,2,3 (deciles becomes Fig 2,
  nested comparison becomes Fig 3); file names are unchanged.
- Supplementary Figs S1-S5 -> order of first citation (weighted S1, binned S2, DCA S3,
  objective AUROC S4, objective calibration S5).
- Reorders SI table blocks and main/SI figure legends accordingly.

Dry run by default; pass --apply to write.
"""
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MS = ROOT / "manuscript"
MANU = MS / "题目4手稿.md"
SI = MS / "Supplementary_Information.md"

REF_MAP = {6:1, 5:2, 9:3, 11:4, 1:5, 2:6, 10:7, 14:8, 15:9, 16:10, 3:11, 12:12, 13:13, 4:14, 7:15, 8:16}
TAB_MAP = {1:1, 12:2, 6:3, 3:4, 5:5, 4:6, 10:7, 11:8, 13:9, 8:10, 9:11, 2:12, 7:13, 14:14}
FIG_MAP_SI = {4:1, 3:2, 1:3, 2:4, 5:5}
MAIN_FIG_SWAP = {1:1, 2:3, 3:2}

def collapse(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j+1] == nums[j] + 1:
            j += 1
        if j - i >= 2: out.append(f"{nums[i]}\u2013{nums[j]}")
        else: out.extend(str(n) for n in nums[i:j+1])
        i = j + 1
    return ",".join(out)

def map_ref_brackets(text):
    def repl(m):
        inner = m.group(1); nums = []
        for part in inner.split(","):
            part = part.strip()
            r = re.match(r"^(\d+)\s*[\u2013-]\s*(\d+)$", part)
            if r: nums += list(range(int(r.group(1)), int(r.group(2)) + 1))
            elif part.isdigit(): nums.append(int(part))
        if not nums: return m.group(0)
        mapped = [REF_MAP.get(n, n) for n in nums]
        return "[" + collapse(mapped) + "]"
    return re.sub(r"\[([0-9,\u2013\- ]+)\]", repl, text)

def map_table_tokens(text):
    # placeholders first, then expand
    text = re.sub(r"S(\d+)(?![\d])", lambda m: f"@@T{int(m.group(1))}@@", text)
    def fin(m): return f"S{TAB_MAP[int(m.group(1))]}"
    return re.sub(r"@@T(\d+)@@", fin, text)

def protect_si_figs(text):
    text = re.sub(r"(Supplementary Figs?\.?\s*)S(\d+)", lambda m: m.group(1) + f"@@F{int(m.group(2))}@@", text)
    return re.sub(r"(Fig\.\s*)S(\d+)", lambda m: m.group(1) + f"@@F{int(m.group(2))}@@", text)

def expand_si_figs(text):
    return re.sub(r"@@F(\d+)@@", lambda m: f"S{FIG_MAP_SI.get(int(m.group(1)), int(m.group(1)))}", text)

def protect_figsr(text):
    return re.sub(r"FigSR(\d+)", lambda m: f"@@FSR{m.group(1)}@@", text)

def restore_figsr(text):
    return re.sub(r"@@FSR(\d+)@@", lambda m: f"FigSR{m.group(1)}", text)

def map_main_fig_labels(text):
    return re.sub(r"(?<!Supplementary )\bFig\.?\s*([123])([A-B])?",
                  lambda m: f"Fig. {MAIN_FIG_SWAP[int(m.group(1))]}" + (m.group(2) or ""), text)

# ---------------- manuscript ----------------
man = MANU.read_text(encoding="utf-8")
head, rest = man.split("## References", 1)
reflist, tail = rest.split("## Figure legends", 1)

# split reference entries
entries = dict(re.findall(r"^(\d+)\.\s+(.*)$", reflist, flags=re.M))
assert len(entries) == 16, len(entries)
new_reflist = "\n".join(f"{new}. {entries[str(old)]}" for old, new in sorted(REF_MAP.items(), key=lambda kv: kv[1]))

# renumber the body (head + figure legends + tables)
head = map_ref_brackets(head)
head = protect_si_figs(head)
head = protect_figsr(head)
head = map_table_tokens(head)
head = expand_si_figs(head)
head = restore_figsr(head)
head = map_main_fig_labels(head)

tail = protect_si_figs(tail)
tail = protect_figsr(tail)
tail = map_table_tokens(tail)
tail = expand_si_figs(tail)
tail = restore_figsr(tail)
tail = map_main_fig_labels(tail)

# ---------------- swap main figure legend blocks 2 and 3 ----------------
lines = tail.splitlines()
blocks = {}
cur = None
buf = []
for ln in lines:
    m = re.match(r"- \*\*Figure ([123]) ", ln)
    if m:
        if cur is not None: blocks[cur] = buf
        cur = int(m.group(1)); buf = [ln]
    elif cur is not None:
        buf.append(ln)
if cur is not None: blocks[cur] = buf
# after mapping, labels are already 1,2,3 but content order is still old 1,2(nested),3(deciles).
# Swap content of new labels 2 and 3.
if 2 in blocks and 3 in blocks:
    b2 = [re.sub(r"\*\*Figure 2", "**Figure 3", x, count=1) for x in blocks[3]]
    b3 = [re.sub(r"\*\*Figure 3", "**Figure 2", x, count=1) for x in blocks[2]]
    blocks[2], blocks[3] = b2, b3
new_tail = "\n".join(blocks[1] + blocks[2] + blocks[3])

new_man = head + "## References\n" + new_reflist + "\n\n## Figure legends\n" + new_tail

# ---------------- SI ----------------
si = SI.read_text(encoding="utf-8")
si = protect_si_figs(si)
si = protect_figsr(si)
si = map_table_tokens(si)
si = expand_si_figs(si)
si = restore_figsr(si)

pre, rest = si.split("## Supplementary tables", 1)
tables_part, legends_part = rest.split("## Supplementary figure legends", 1)
# split blocks
parts = re.split(r"(?=^### Supplementary Table S\d+\.)", tables_part, flags=re.M)
prefix = parts[0]
blocks = parts[1:]
def num(b): return int(re.search(r"### Supplementary Table S(\d+)\.", b).group(1))
blocks.sort(key=num)
new_tables = prefix.rstrip() + "\n\n" + "\n\n".join(b.strip() for b in blocks) + "\n\n"

# SI legends: split and reorder by new label
leg_lines = legends_part.strip().splitlines()
leg_header = leg_lines[0] if leg_lines and leg_lines[0].startswith("##") else "## Supplementary figure legends"
leg_blocks = {}
cur = None; buf = []
_start = 1 if leg_lines and leg_lines[0].startswith("##") else 0
for ln in leg_lines[_start:]:
    m = re.match(r"- \*\*Supplementary Fig\. S(\d+)", ln)
    if m:
        if cur is not None: leg_blocks[cur] = buf
        cur = int(m.group(1)); buf = [ln]
    elif cur is not None:
        buf.append(ln)
if cur is not None: leg_blocks[cur] = buf
new_legends = leg_header + "\n" + "\n".join("\n".join(leg_blocks[k]) for k in sorted(leg_blocks))
new_si = pre + "## Supplementary tables" + new_tables + new_legends + "\n"

if "--apply" in sys.argv:
    MANU.write_text(new_man, encoding="utf-8")
    SI.write_text(new_si, encoding="utf-8")

print("refs new order:", [new for old,new in sorted(REF_MAP.items(), key=lambda kv: kv[1])])
print("table order in SI:", [num(b) for b in blocks])
print("figure legend order:", sorted(leg_blocks))
print([ln[:70] for ln in new_si.splitlines() if ln.startswith("- **Supplementary Fig")])
print("DCA lines:",[ln[:100] for ln in new_si.splitlines() if "DCA" in ln])
print("applied:", "--apply" in sys.argv)
print("--- first 600 chars of new refs ---")
print(new_reflist[:600])