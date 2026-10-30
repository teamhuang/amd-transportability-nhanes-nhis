import pandas as pd, numpy as np
def probe(y):
    df = pd.read_csv(rf"H:\眼科论文\题目4\data\raw\nhis\e{y}\samadult.csv")
    df.columns=[c.upper() for c in df.columns]
    a = df[df['AGE_P'].between(40,999)]
    print(f"\n===== NHIS {y}: rows={len(df)}, 40+={len(a)} =====")
    for v in ['SMKSTAT2','CIGAREV2','CIGCUR2']:
        if v in df.columns: print(v, dict(df[v].value_counts(dropna=False).sort_index().head(6)))
    def yes(var, eq):
        s=pd.to_numeric(a[var],errors='coerce'); return (s==eq).sum(), len(a)
    for var,code in [('VIM_MDEV',1),('HYPEV',1),('CHLEV',1),('DIBEV1',1),('VIMGLASS',1)]:
        if var in a:
            n,tt=yes(var,code); print(f"40+ {var}=={code}: {n} ({n/1e4 if False else 0:.0%}) prev={n/len(a):.4f}")
    # smoking status among 40+
    if 'SMKSTAT2' in a: print("40+ SMKSTAT2:", dict(pd.to_numeric(a['SMKSTAT2'],errors='coerce').value_counts(dropna=False).sort_index()))
    # weight availability
    for w in ['WTFA_SA','WTIA_SA']:
        if w in a:
            ww=pd.to_numeric(a[w],errors='coerce'); print(w,"nonnull",int(ww.notna().sum()),"median",round(float(ww.median()),1) if ww.notna().any() else None)
    print("PSU/STRAT present:", [c for c in a.columns if c in ('PPSU','PSTRAT','REGION')])
for y in [2016,2017]:
    probe(y)
