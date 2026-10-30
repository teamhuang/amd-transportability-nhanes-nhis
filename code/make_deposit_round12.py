"""Round-12 P0-2: assemble a submission-ready code/data deposit bundle.

Creates a staging tree + zip that the author can upload to GitHub and archive on Zenodo.
It deliberately does NOT invent a DOI, licence or repository URL: those are author-only
decisions, flagged in README_DEPOSIT.md.

Layout:
  code/            final analysis scripts (plus development probes in code/dev/)
  results/         every JSON + coefficient table + derived CSVs emitted by the analyses
  derived_data/    the two analytic datasets used for all reported numbers
  documentation/   harmonisation dictionary, endpoint audit log, TRIPOD checklist, design freeze
  figures/         manuscript + supplementary figures (PDF/PNG)
  environment/     requirements.txt + environment.txt
  MANIFEST.json    SHA-256 and row counts for every deposited file
"""
import hashlib, json, platform, shutil, sys, zipfile
from importlib import metadata
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SCR = ROOT / "analysis" / "scripts"
OUT = ROOT / "analysis" / "output"
MS = ROOT / "manuscript"
FIG = ROOT / "figures"
STAGE = Path(r"C:\Users\Admin\.qoderworkcn\workspace\mvs2di6b7792q6ts\deposit_staging")
ZIP = Path(r"C:\Users\Admin\.qoderworkcn\workspace\mvs2di6b7792q6ts\outputs\topic4_amd_transport_deposit_v1.zip")
ZIP.parent.mkdir(parents=True, exist_ok=True)

if STAGE.exists():
    shutil.rmtree(STAGE)          # QoderWork scratch staging tree, rebuilt from scratch each run
for sub in ["code", "code/dev", "results", "derived_data", "documentation", "figures", "environment"]:
    (STAGE / sub).mkdir(parents=True, exist_ok=True)

manifest = {"created": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "files": []}

def add(src: Path, rel: str):
    dst = STAGE / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    h = hashlib.sha256(dst.read_bytes()).hexdigest()
    entry = {"path": rel, "bytes": dst.stat().st_size, "sha256": h}
    if dst.suffix.lower() == ".csv":
        try:
            entry["rows"] = int(len(pd.read_csv(dst, low_memory=False)))
        except Exception:
            entry["rows"] = None
    manifest["files"].append(entry)

# code
primary, dev = [], []
for p in sorted(SCR.glob("*.py")):
    if p.name.startswith("_"):
        dev.append(p)
    else:
        primary.append(p)
for p in primary:
    add(p, f"code/{p.name}")
for p in dev:
    add(p, f"code/dev/{p.name}")

# results — outputs of the superseded derive_and_validate.py fit go to results/legacy/
LEGACY = {"model_coefficients.csv", "external_lp.csv", "external_metrics_full.json",
          "external_metrics_flat.csv", "derivation_pred.csv", "age_dose.json",
          "findings.json", "internal_validation.json", "score_points.json",
          "sensitivity_objective.json"}
for p in sorted(OUT.glob("*")):
    if p.is_file() and p.suffix.lower() in {".json", ".csv", ".txt", ".md"}:
        sub = "results/legacy" if p.name in LEGACY else "results"
        add(p, f"{sub}/{p.name}")

legacy_note = """# Legacy outputs — do not use to reproduce manuscript numbers

These files were produced by `code/derive_and_validate.py`, an **earlier fit that used a
restricted-cubic-spline age term and an earlier coefficient set**. They are retained only for audit
history. The authoritative, final (linear-age) outputs are the `sr_*.json` files in the parent
directory, and the authoritative coefficients are in `model_equation.json`
(`coefficients_full_precision`).

Concrete examples of the discrepancy:
- `model_coefficients.csv` lists spline terms `cr(age, df=3)[0..2]`; the published model has a single
  linear `age` coefficient (0.09447).
- `external_metrics_full.json` / `external_lp.csv` give a pooled NHIS AUROC of 0.7748, mean predicted
  4.43% and slope 0.848, whereas the published transport gives 0.7737, 4.36% and 0.8723
  (`sr_external_age.json`).

If these files disagree with the manuscript, the `sr_*` files and `model_equation.json` are correct.
"""
(STAGE / "results" / "legacy").mkdir(parents=True, exist_ok=True)
(STAGE / "results" / "legacy" / "README_LEGACY.md").write_text(legacy_note, encoding="utf-8")
manifest["files"].append({"path": "results/legacy/README_LEGACY.md",
                          "bytes": (STAGE / "results" / "legacy" / "README_LEGACY.md").stat().st_size,
                          "sha256": hashlib.sha256((STAGE / "results" / "legacy" / "README_LEGACY.md").read_bytes()).hexdigest()})

# derived data
for name in ["nhanes_derivation.csv", "nhis_external.csv"]:
    add(OUT / name, f"derived_data/{name}")

# documentation
for name in ["harmonization_dictionary.md", "endpoint_code_audit.md",
             "TRIPOD_checklist.md", "Submission_Readiness_Checklist.md",
             "Supplementary_Information.md", "reference_verification.md", "Round16_Change_Log.md"]:
    add(MS / name, f"documentation/{name}")
add(ROOT / "admission_record.md", "documentation/admission_record.md")

# figures
for p in sorted(FIG.glob("*")):
    if p.is_file() and p.suffix.lower() in {".pdf", ".png"}:
        add(p, f"figures/{p.name}")

# environment
vers = {}
for pkg in ["pandas", "numpy", "statsmodels", "scikit-learn", "scipy", "matplotlib",
            "pyreadstat", "patsy", "pypdf"]:
    try:
        vers[pkg] = metadata.version(pkg)
    except Exception:
        vers[pkg] = "not installed"
req = "\n".join(f"{k}=={v}" for k, v in vers.items()) + "\n"
(STAGE / "environment" / "requirements.txt").write_text(req, encoding="utf-8")
(STAGE / "environment" / "environment.txt").write_text(
    f"Python {sys.version}\nPlatform: {platform.platform()}\n\n" +
    "\n".join(f"{k}=={v}" for k, v in vers.items()) + "\n", encoding="utf-8")
for name in ["requirements.txt", "environment.txt"]:
    add(STAGE / "environment" / name, f"environment/{name}")

readme = f"""# Deposit — no-imaging model for self-reported AMD: derivation (NHANES 2005–08) and external validation (NHIS 2016/2017)

Assembled {manifest['created']} from the final analysis state. Every reported number in the
manuscript is emitted by the scripts in `code/` and stored in `results/`.

## Contents
- `code/` — analysis scripts. Primary chain: `build_nhanes.py` → `build_nhis.py` → `harden.py` →
  `harden2.py` → `sr_derivation.py` → `sr_weighted_calibration.py` → `sr_survey_inference.py` →
  `sr_increment_robust.py` → `sr_missing_mi2.py` → `sr_mi_stability.py` → `sr_lipid_sensitivity.py` →
  `sr_final_tables_figs.py` → `make_tables_sr.py` → `make_s6_table.py` → `make_round4_tables_figs.py`.
  Verification/audit: `verify_round2.py`, `verify_round3.py`, `verify_round4.py`, `audit_round2.py`,
  `audit_round3.py`, `audit_round5.py`, `audit_round7.py`, `audit_survey_round8.py`,
  `audit_codebook_round12.py`, `regen_tables_round12.py`, `split_main_si_round15.py`. `code/dev/` holds
  exploratory probes and superseded scripts (`final_equation.py`, `make_tables.py`, `make_tables2.py`)
  kept for transparency.
- **Submission structure:** the main manuscript carries 3 figures and 5 tables (8 display items, the
  Scientific Reports limit); Supplementary Tables S1–S14 and Supplementary Figs S1–S5 are in
  `documentation/Supplementary_Information.md`.
- `results/` — all machine-readable outputs (JSON/CSV) and run logs, including the authoritative
  full-precision coefficients (`model_equation.json`, key `coefficients_full_precision`), the codebook
  audit (`sr_codebook_audit.json`) and the table-regeneration check (`sr_table_regen.json`).
  **`results/legacy/` holds outputs of the earlier spline-model fit; see its README and do not use them
  to reproduce the manuscript.**
- `derived_data/` — the two analytic datasets used for all reported numbers
  (`nhanes_derivation.csv`, `nhis_external.csv`). Raw NCHS transport files are **not** redistributed;
  download them from the CDC/NCHS portals (see manuscript Data availability).
- `documentation/` — harmonisation dictionary, endpoint-code audit log, TRIPOD+AI checklist,
  admission record (design freeze + decision log), submission-readiness checklist.
- `figures/` — manuscript and supplementary figures (PDF + PNG).
- `environment/` — exact package versions.

## How to reproduce
1. Install `environment/requirements.txt` into Python {manifest['python']}.
2. Download the NHANES 2005–06/2007–08 and NHIS 2016/2017 public-use files from CDC/NCHS and place them
   where `build_nhanes.py` / `build_nhis.py` expect (paths are defined at the top of each script).
3. Run the primary chain in the order above; then `regen_tables_round12.py` to prove the display tables
   match the manuscript line-by-line.

## Author decisions still required before upload
- Choose a code licence (e.g. MIT / BSD-3) and add `LICENSE`.
- Create the GitHub repository and the Zenodo **version** DOI, then insert the DOI in the manuscript's
  Code availability section and in this README.
- Confirm that the derived datasets may be redistributed (they contain only de-identified public-use
  variables, but this is the authors' decision).

## Integrity
`MANIFEST.json` lists the SHA-256 and byte size of every file, plus row counts for CSV files.
"""
(STAGE / "README_DEPOSIT.md").write_text(readme, encoding="utf-8")
manifest["files"].append({
    "path": "README_DEPOSIT.md",
    "bytes": (STAGE / "README_DEPOSIT.md").stat().st_size,
    "sha256": hashlib.sha256((STAGE / "README_DEPOSIT.md").read_bytes()).hexdigest()})
(STAGE / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

if ZIP.exists():
    ZIP.unlink()          # generated artifact inside the QoderWork workspace, regenerated on each run
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    z.write(STAGE / "MANIFEST.json", "MANIFEST.json")
    z.write(STAGE / "README_DEPOSIT.md", "README_DEPOSIT.md")
    for p in sorted(STAGE.rglob("*")):
        if p.is_file() and p.name not in {"MANIFEST.json", "README_DEPOSIT.md"}:
            z.write(p, str(p.relative_to(STAGE)))

summary = {"zip": str(ZIP), "zip_mb": round(ZIP.stat().st_size / 1e6, 1),
           "n_files": len(manifest["files"]), "python": manifest["python"],
           "packages": vers,
           "rows": {f["path"]: f["rows"] for f in manifest["files"] if f.get("rows") is not None}}
(STAGE / "DEPOSIT_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
