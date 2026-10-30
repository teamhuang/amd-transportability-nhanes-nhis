import pypdf, re, sys
path = sys.argv[1]
toks = ["VIM_MDEV","HYPEV","CHLEV","DIBEV1","BMI","SMOK","CIG","AGE_P","SEX","REGION","PSU","STRAT","PERW","WTMEC","PERINT","WTA","WTFA","SMOKSTAT","SMOKED","CIGAREV","CIGCUR","HEIGHT","WEIGHT"]
r = pypdf.PdfReader(path)
full = []
for p in r.pages:
    t = p.extract_text() or ""
    full.append(t)
txt = "\n".join(full)
lines = txt.splitlines()
print("TOTAL LINES", len(lines), "PAGES", len(r.pages))
for tk in toks:
    hits = [ln.strip() for ln in lines if re.search(r'\b'+re.escape(tk)+r'\b', ln, re.I)]
    print(f"\n#### {tk} ({len(hits)}) ####")
    for h in hits[:6]:
        print("  ", h[:140])
