import pyreadstat, os, sys
RAW = r"H:\眼科论文\眼科2\data\raw\nhanes"
def labels(f):
    df, meta = pyreadstat.read_xport(os.path.join(RAW, f))
    print("===", f, "rows", len(df))
    m = meta.column_names_to_labels
    for v in meta.column_names:
        print(v, "|", (m.get(v) or "")[:90])
for f in ["VIQ_D.xpt"]:
    labels(f)
print("\n\n#### OPXRET key grading vars ####")
df, meta = pyreadstat.read_xport(os.path.join(RAW, "OPXRET_E.xpt"))
m2 = meta.column_names_to_labels
for v in meta.column_names:
    if v.upper().startswith("OPDD") or "ARM" in v.upper():
        print(v, "|", (m2.get(v) or "")[:80])
