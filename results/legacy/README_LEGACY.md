# Legacy outputs — do not use to reproduce manuscript numbers

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
