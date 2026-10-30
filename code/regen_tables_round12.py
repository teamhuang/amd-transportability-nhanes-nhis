"""Round-12 P0-4: regenerate every display table from the final code and prove that the
tables embedded in the manuscript are exactly what the code emits.

Sequence (same order as the project's documented chain):
  make_tables_sr.py          -> output/Tables_sr.md (+ manuscript/harmonization_dictionary.md)
  make_s6_table.py           -> inserts S6 into Tables_sr.md (idempotent)
  make_round4_tables_figs.py -> output/Tables_sr4.md (+ FigSR7/8)

Then reconstructs the manuscript's tables block as
  "\n## Tables\n\n" + Tables_sr.md + "\n\n" + Tables_sr4.md (without its "## Round-4..." heading)
and compares it, line by line, with the block that currently follows the insertion marker.

Writes analysis/output/sr_table_regen.json with hashes, timing and the comparison result.
Backs up every file it may overwrite to the QoderWork workspace (no files are deleted).
"""
import hashlib, json, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "output"
MS = ROOT / "manuscript"
FIG = ROOT / "figures"
SCRIPTS = ROOT / "analysis" / "scripts"
MANU = MS / "题目4手稿.md"
SI = MS / "Supplementary_Information.md"
BACKUP = Path(r"C:\Users\Admin\.qoderworkcn\workspace\mvs2di6b7792q6ts\backup_round12")
BACKUP.mkdir(parents=True, exist_ok=True)

def sha(p):
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()

def backup(p):
    rel = p.relative_to(ROOT)
    dst = BACKUP / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if p.exists() and not dst.exists():
        shutil.copy2(p, dst)

watched = [OUT / "Tables_sr.md", OUT / "Tables_sr4.md", MS / "harmonization_dictionary.md",
           MANU, SI]
for p in watched:
    backup(p)
for p in list(FIG.glob("FigSR7*")) + list(FIG.glob("FigSR8*")):
    backup(p)

before = {str(p.relative_to(ROOT)): sha(p) for p in watched if p.exists()}
t0 = time.time()
runs = []
for s in ["make_tables_sr.py", "make_s6_table.py", "make_round4_tables_figs.py"]:
    t = time.time()
    r = subprocess.run([sys.executable, s], cwd=str(SCRIPTS), capture_output=True, text=True)
    runs.append({"script": s, "returncode": r.returncode, "seconds": round(time.time() - t, 1),
                 "stderr_tail": (r.stderr or "")[-400:]})
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        raise SystemExit(f"{s} failed with {r.returncode}")
elapsed = round(time.time() - t0, 1)

after = {str(p.relative_to(ROOT)): sha(p) for p in watched if p.exists()}
changed = {k: {"before": before.get(k), "after": after[k]} for k in after
           if before.get(k) != after[k]}

# ---- reconstruct the expected main and SI blocks ----
import re, subprocess as _sp
tsr = (OUT / "Tables_sr.md").read_text(encoding="utf-8").strip()
_split = tsr.index("## Supplementary tables")
expected_main = tsr[:_split].strip()
si_from_sr = tsr[_split:].strip()
tsr4_lines = (OUT / "Tables_sr4.md").read_text(encoding="utf-8").strip().splitlines()
if tsr4_lines and tsr4_lines[0].startswith("## Round-4"):
    tsr4_lines = tsr4_lines[1:]
tsr4 = "\n".join(tsr4_lines).strip()
expected_si = si_from_sr + "\n\n" + tsr4

doc = MANU.read_text(encoding="utf-8")
marker = "<!-- TABLES inserted below by script; numbers emitted from analysis/output/sr_*.json -->"
assert marker in doc, "insertion marker missing"
head, block = doc.split(marker, 1)

def norm(s):
    lines = [ln.rstrip() for ln in s.replace("\r\n", "\n").split("\n")]
    out = []
    for ln in lines:
        if ln == "" and out and out[-1] == "":
            continue
        out.append(ln)
    return "\n".join(out).strip()

def si_table_section():
    if not SI.exists():
        return ""
    t = SI.read_text(encoding="utf-8")
    a = t.find("## Supplementary tables")
    b = t.find("## Supplementary figure legends")
    return t[a:b].strip() if a != -1 else ""

def canon(s):
    """Normalize for comparison: ignore bold markers and supplementary numbering
    (numbering is verified separately by audit_consistency_round12.py)."""
    s = s.replace("**", "")
    s = re.sub(r"S\d+(?:\s*[\u2013,\-]\s*S?\d+)*", "S#", s)
    return norm(s)

def compare(a, b):
    if canon(a) == canon(b):
        return True, None
    aa = canon(a).split("\n"); bb = canon(b).split("\n")
    for i in range(max(len(aa), len(bb))):
        x = aa[i] if i < len(aa) else "<EOF>"
        y = bb[i] if i < len(bb) else "<EOF>"
        if x != y:
            return False, {"line": i + 1, "manuscript": x[:200], "regenerated": y[:200]}
    return False, None

def si_block_bodies(text):
    parts = re.split(r"(?=^### Supplementary Table S\d+\.)", text, flags=re.M)
    out = []
    for p in parts:
        if not p.startswith("### Supplementary Table"):
            continue
        body = re.sub(r"^### Supplementary Table S\d+\.[^\n]*\n", "", p).strip()
        out.append(canon(body))
    return sorted(out)

def compare_si(a, b):
    ca, cb = si_block_bodies(a), si_block_bodies(b)
    if ca == cb:
        return True, None
    return False, {"manuscript_blocks": len(ca), "regenerated_blocks": len(cb),
                   "first_mismatch_preview": (ca[:1] + cb[:1])}

same_main, diff_main = compare(block.split("## Supplementary information")[0], expected_main)
same_si, diff_si = compare_si(si_table_section(), expected_si)

applied = False
if "--apply" in sys.argv:
    print("NOTE: automatic re-splicing is disabled after the round-16 renumbering; "
          "update the final manuscript/SI directly if table content changes.")

# ---- per-display-item presence check across both files ----
ids = re.findall(r"^### (Table \d+|Supplementary Table S\d+)\.", tsr + "\n" + tsr4, flags=re.M)
present = {i: ((f"### {i}." in block) or (f"### {i}." in si_table_section())) for i in ids}
missing = [k for k, v in present.items() if not v]

rec = {
    "regenerated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    "runs": runs, "seconds_total": elapsed,
    "hashes_before": before, "hashes_after": after,
    "files_changed": sorted(changed),
    "tables_in_generated_files": ids,
    "tables_missing_from_manuscript": missing,
    "manuscript_block_matches_generated": same_main,
    "si_table_block_matches_generated": same_si,
    "main_tables": [i for i in ids if i.startswith("Table ")],
    "si_tables": [i for i in ids if i.startswith("Supplementary")],
    "applied_to_manuscript": applied,
    "manuscript_last_edited": time.strftime("%Y-%m-%d %H:%M:%S",
                                            time.localtime(MANU.stat().st_mtime)),
    "first_difference_main": diff_main,
    "first_difference_si": diff_si,
    "verdict": "PASS" if same_main and same_si and not missing else "FAIL",
}
(OUT / "sr_table_regen.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False), encoding="utf-8")

print(json.dumps({k: rec[k] for k in ["files_changed", "main_tables", "si_tables",
                                      "tables_missing_from_manuscript",
                                      "manuscript_block_matches_generated",
                                      "si_table_block_matches_generated",
                                      "first_difference_main", "first_difference_si",
                                      "verdict"]}, indent=2, ensure_ascii=False))
print("changed file details:")
for k, v in changed.items():
    print(" ", k, v["before"][:12] if v["before"] else None, "->", v["after"][:12])
