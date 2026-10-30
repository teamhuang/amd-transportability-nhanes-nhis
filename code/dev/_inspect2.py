import pyreadstat, os
import numpy as np
RAW = r"H:\眼科论文\眼科2\data\raw\nhanes"
def rd(f):
    df, meta = pyreadstat.read_xport(os.path.join(RAW, f))
    df.columns=[c.upper() for c in df.columns]
    return df
bpx = rd("BPX_E.xpt")
print("BPX cols:", [c for c in bpx.columns if c.startswith("BPX")])
opx = rd("OPXRET_E.xpt")
for v in ["OPDUARMA","OPDUARM","OPDULARM","OPDUARMX","OPDDARM","OPDDLARM","OPDDEXU","OPDDAGA"]:
    if v in opx.columns:
        print(v, dict(opx[v].value_counts(dropna=False).sort_index()))
print("TRIGLY cols:", [c for c in rd("TRIGLY_E.xpt").columns])
print("TCHOL cols:", [c for c in rd("TCHOL_E.xpt").columns])
print("HDL cols:", [c for c in rd("HDL_E.xpt").columns])
print("SMQ cols:", [c for c in rd("SMQ_E.xpt").columns if c.startswith("SMQ")][:8])
