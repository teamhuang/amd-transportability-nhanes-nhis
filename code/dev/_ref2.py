import sys; sys.path.insert(0,'analysis/scripts')
from verify_refs import pubmed
qs={
 "tripod2015":"Collins Reitsma Altman TRIPOD 2015 BMJ reporting clinical prediction models",
 "amd_prevalence_us":"Friedman Prevalence of age-related macular degeneration United States Archives Ophthalmology 2004",
 "eye_case_control":"Risk factors for age-related macular degeneration Eye Disease Case-Control Study Group 1992",
 "blue_mountains":"Wang Mitchell Cumunan age-related macular degeneration Blue Mountains Eye Study prevalence",
 "amd_lancet":"Age-related macular degeneration Lancet review Bird",
 "nhis_health_status":"National Health Interview Survey vision eye health prevalence data",
 "areds_2":"Age-Related Eye Disease Study design baseline AREDS Report No 1",
}
for k,q in qs.items():
    print('\n#####',k)
    for r in pubmed(q)[:2]:
        print(f"[PMID {r['pmid']}] {r['authors']} ({r['year']}) {r['journal']}. {r['title']}")
