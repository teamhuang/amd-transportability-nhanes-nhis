import sys; sys.path.insert(0,'analysis/scripts')
from verify_refs import pubmed
qs={
 "global_burden":"Global prevalence age-related macular degeneration and disease burden projection 2020 2040",
 "friedman":"Friedman Prevalence age-related macular degeneration United States",
 "prediction_rules":"methodologic rules for clinical prediction models Steyerberg Tripepi",
 "nhanes_vision":"National Health and Nutrition Examination Survey eye health vision impairment prevalence",
 "nhis_amd":"National Health Interview Survey macular degeneration prevalence United States",
 "amd_polygenic":"polygenic risk score prediction age-related macular degeneration",
 "amd_ml":"machine learning deep learning prediction age-related macular degeneration fundus images external validation",
}
for k,q in qs.items():
    print('\n#####',k)
    for r in pubmed(q)[:2]:
        print(f"[PMID {r['pmid']}] {r['authors']} ({r['year']}) {r['journal']}. {r['title']}")
