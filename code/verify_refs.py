import urllib.request, urllib.parse, json, time, sys
def get(url):
    req=urllib.request.Request(url,headers={"User-Agent":"ref-check/1.0"})
    return urllib.request.urlopen(req,timeout=30).read().decode("utf-8","ignore")
def pubmed(term):
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=3&term="+urllib.parse.quote(term)
    ids=json.loads(get(u))['esearchresult']['idlist']
    if not ids: return []
    s=json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="+",".join(ids)))
    out=[]
    for i in ids:
        d=s['result'][i]
        auth=", ".join([a['name'] for a in d.get('authors',[])][:6])
        out.append(dict(pmid=i,title=d.get('title',''),journal=d.get('source',''),year=d.get('pubdate','')[:4],authors=auth))
    return out
def crossref(q):
    u="https://api.crossref.org/works?rows=1&query.bibliographic="+urllib.parse.quote(q)
    try:
        j=json.loads(get(u))['message']['items'][0]
        return dict(doi=j.get('DOI'),title=(j.get('title') or [''])[0],year=(j.get('issued',{}).get('date-parts',[[None]])[0][0]),journal=(j.get('container-title') or [''])[0])
    except Exception as e:
        return dict(doi=None,err=str(e))

queries={
 "tripod":"TRIPOD transparent reporting multivariable prediction model Collins 2015 BMJ",
 "tripod_ai":"TRIPOD-AI artificial intelligence prediction models Collins Moons 2025",
 "calib":"Van Calster calibration Achilles heel predictive models medical decision making",
 "dca":"Vickers Elkin decision curve analysis medical decision making 2006",
 "amd_epi":"Wong age-related macular degeneration Lancet review",
 "amd_risk_areds":"risk factors age-related macular degeneration AREDS eye disease",
 "blue_mountain":"Blue Mountains Study age-related macular degeneration risk Wang Mitchell",
 "nhis":"National Health Interview Survey description data",
 "nhanes":"National Health and Nutrition Examination Survey plan analysis",
 "steyerberg":None,
}
for k,q in queries.items():
    if q is None: continue
    print("\n#####",k)
    try:
        pm=pubmed(q)
        for r in pm[:2]:
            cr=crossref(r['title'][:80])
            print(f"[PMID {r['pmid']}] {r['authors']} ({r['year']}) {r['journal']}. {r['title']} | DOI:{cr.get('doi')}")
    except Exception as e:
        print("err",e)
    time.sleep(0.4)
