"""Round-12/15 consistency audit: main manuscript + Supplementary Information.

Checks figure/table cross-references (including Supplementary Fig. S1-S5), the main display-item
count, minus-sign convention, decimal precision, endpoint terminology and revision language.
Reports only; no file is modified.
"""
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MS = ROOT / "manuscript"
FIG = ROOT / "figures"
OUT = ROOT / "analysis" / "output"
man = (MS / "题目4手稿.md").read_text(encoding="utf-8")
si_path = MS / "Supplementary_Information.md"
si = si_path.read_text(encoding="utf-8") if si_path.exists() else ""
combined = man + "\n" + si
MARKER = "<!-- TABLES inserted below by script; numbers emitted from analysis/output/sr_*.json -->"
res = {"checks": [], "findings": []}

def add(name, ok, detail):
    res["checks"].append({"check": name, "pass": bool(ok), "detail": detail})

# ---- 1. main figures ----
main_figs = sorted({int(m.group(1)) for m in re.finditer(r"(?<!Supplementary )\bFig(?:ure)?s?\s?(\d+)\b", man)})
missing_main = [n for n in main_figs if not list(FIG.glob(f"Fig{n}_*.pdf"))]
add("every main figure reference resolves to a file", not missing_main, str(missing_main))
res["main_figures"] = main_figs

# ---- 1b. supplementary figures ----
si_figs = sorted({int(m.group(1)) for m in re.finditer(r"Supplementary Fig(?:ure)?s?\.?\s?S(\d+)", combined)})
missing_si = [n for n in si_figs if not list(FIG.glob(f"FigSR{n}_*.pdf"))]
legend_missing = [n for n in si_figs if f"Supplementary Fig. S{n}" not in si]
add("every supplementary figure reference resolves to a file", not missing_si, str(missing_si))
add("every supplementary figure has an SI legend", not legend_missing, str(legend_missing))
res["si_figures"] = si_figs

# ---- 1c. display-item count (main text) ----
n_main_tables = len(re.findall(r"^### Table \d+\.", man, flags=re.M))
add("main display items <= 8", len(main_figs) + n_main_tables <= 8,
    f"{len(main_figs)} figures + {n_main_tables} tables = {len(main_figs)+n_main_tables}")

# ---- 1e. numbering by order of first appearance ----
body = man.split("## References")[0]

def expand_s_list(chunk):
    nums = []
    for a, b in re.findall(r"S(\d+)(?:\s*[\u2013-]\s*S?(\d+))?", chunk):
        a = int(a); b = int(b) if b else a
        nums += list(range(a, b + 1))
    return nums

def first_order(seq):
    out = []
    for n in seq:
        if n not in out: out.append(n)
    return out

# references
ref_seq = []
for m in re.finditer(r"\[([0-9,\u2013\- ]+)\]", body):
    for part in m.group(1).split(","):
        part = part.strip()
        r = re.match(r"^(\d+)\s*[\u2013-]\s*(\d+)$", part)
        if r: ref_seq += list(range(int(r.group(1)), int(r.group(2)) + 1))
        elif part.isdigit(): ref_seq.append(int(part))
ref_order = first_order(ref_seq)
ref_ok = ref_order == list(range(1, 17))
add("references numbered by first appearance", ref_ok, str(ref_order))

# supplementary tables
stab_seq = []
for m in re.finditer(r"Supplementary Tables?((?:\s*S\d+(?:\s*[\u2013,-]\s*S?\d+)?)+)", body):
    stab_seq += expand_s_list(m.group(1))
stab_order = first_order(stab_seq)
stab_ok = stab_order == list(range(1, 15))
add("supplementary tables numbered by first appearance", stab_ok, str(stab_order))

# main figures
main_seq = [int(m.group(1)) for m in re.finditer(r"(?<!Supplementary )\bFig\.?\s*(\d+)", body)]
main_order = first_order(main_seq)
main_ok = main_order == [1, 2, 3]
add("main figures numbered by first appearance", main_ok, str(main_order))

# supplementary figures
sifig_seq = []
for m in re.finditer(r"Supplementary Figs?\.?\s*S(\d+(?:\s*[,\u2013-]\s*S?\d+)*)", body):
    for part in re.split(r"[,\u2013-]", m.group(1)):
        part = part.strip().lstrip("S")
        if part.isdigit(): sifig_seq.append(int(part))
sifig_order = first_order(sifig_seq)
sifig_ok = sifig_order == [1, 2, 3, 4, 5]
add("supplementary figures numbered by first appearance", sifig_ok, str(sifig_order))

# every display item cited
tables_cited = set(stab_order)
missing_table_cites = [i for i in range(1, 15) if i not in tables_cited]
missing_main_fig_cites = [i for i in range(1, 4) if i not in main_order]
missing_si_fig_cites = [i for i in range(1, 6) if i not in sifig_order]
add("every supplementary table is cited", not missing_table_cites, str(missing_table_cites))
add("every main figure is cited", not missing_main_fig_cites, str(missing_main_fig_cites))
add("every supplementary figure is cited", not missing_si_fig_cites, str(missing_si_fig_cites))

# ---- 1f. bold only on figure titles ----
bad_bold = []
for ln in (man + "\n" + si).splitlines():
    if "**" in ln and not (ln.startswith("- **Figure ") or ln.startswith("- **Supplementary Fig.")):
        bad_bold.append(ln[:120])
add("bold used only in figure titles", not bad_bold, str(bad_bold[:5]))

# ---- 1d. figure file integrity ----
bad_files = []
inventory = []
for pdf in sorted(FIG.glob("*.pdf")):
    size = pdf.stat().st_size; head = pdf.read_bytes()[:5]
    ok = head == b"%PDF-" and size > 10_000
    inventory.append({"file": pdf.name, "bytes": size, "ok": ok})
    if not ok: bad_files.append(pdf.name)
for png in sorted(FIG.glob("*.png")):
    size = png.stat().st_size; head = png.read_bytes()[:8]
    ok = head == b"\x89PNG\r\n\x1a\n" and size > 10_000
    inventory.append({"file": png.name, "bytes": size, "ok": ok})
    if not ok: bad_files.append(png.name)
add("every figure file is a valid, non-trivial PDF/PNG", not bad_files, str(bad_files))
res["figure_files"] = inventory

# ---- 2. tables ----
headings = set(re.findall(r"^### (Table \d+|Supplementary Table S\d+)\.", combined, flags=re.M))
refs = set()
for m in re.finditer(r"\bSupplementary Table S(\d+)", combined):
    refs.add(f"Supplementary Table S{m.group(1)}")
for m in re.finditer(r"(?<!Supplementary )\bTable\s+(\d+)", combined):
    refs.add(f"Table {m.group(1)}")
missing_tabs = sorted(r for r in refs if r not in headings)
add("every table reference resolves to a heading", not missing_tabs, str(missing_tabs))
res["table_headings"] = sorted(headings, key=lambda x: (x.startswith("Supplementary"), x))
res["table_refs"] = sorted(refs)

# ---- 3. minus signs: narrative vs tables ----
narr, tabs = man.split(MARKER, 1) if MARKER in man else (man, "")
def minus_counts(s):
    return {"u2212": len(re.findall(r"\u2212\s?\d", s)),
            "ascii": len(re.findall(r"(?<![\w.])-\s?\d+\.\d", s))}
mn, mt = minus_counts(narr), minus_counts(tabs + "\n" + si)
add("narrative uses one minus-sign convention", min(mn["u2212"], mn["ascii"]) == 0, str(mn))
add("tables use one minus-sign convention", min(mt["u2212"], mt["ascii"]) == 0, str(mt))
res["minus_sign"] = {"narrative": mn, "tables": mt}

# ---- 4. decimal precision in the narrative ----
body = man.split("## References")[0]
threes = len(re.findall(r"\b0\.\d{3}\b", body))
fours = len(re.findall(r"\b0\.\d{4}\b", body))
res["decimal_precision"] = {"3dp": threes, "4dp": fours}
add("decimal precision is deliberate", fours <= threes, f"3 dp={threes}, 4 dp={fours}")

# ---- 5. endpoint terminology ----
terms = {t: len(re.findall(t, combined, flags=re.I)) for t in
         ["any-ARM", "any ARM", "late AMD", "late ARM", "early ARM", "photographic"]}
res["terminology_counts"] = terms
bad = [m.group(0).strip()[:160] for m in re.finditer(r"[^.]*\blate AMD\b[^.]*\.", combined)
       if "OPDULARM" not in m.group(0) and "fundus" not in m.group(0) and "photograph" not in m.group(0)]
add("'late AMD' always tied to the photographic OPDULARM definition", not bad, str(bad[:5]))

# ---- 6. revision language ----
banned = ["reviewer flagged", "reviewer request", "as requested by the reviewer", "v8-SR", "v9-SR",
          "the reviewer", "response to review", "(words)", "word count", "load-bearing",
          "heterogeneity survived", "The honest conclusion", "extreme bounds", "was protective"]
hits = [b for b in banned if b.lower() in combined.lower()]
add("no revision-process language", not hits, str(hits))

res["verdict"] = "PASS" if all(c["pass"] for c in res["checks"]) else "REVIEW"
(OUT / "sr_consistency_round12.json").write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(res, indent=2, ensure_ascii=False))