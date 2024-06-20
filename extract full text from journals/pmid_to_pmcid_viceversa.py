import requests

def pmc_to_pmid(pmcid_list):
    pmcids = ','.join(pmcid_list)
    url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&ids={pmcids}"
    response = requests.get(url)
    
    if response.status_code == 200:
        results = response.json()
        pmid_dict = {}
        for record in results['records']:
            if 'pmid' in record:
                pmid_dict[record['pmcid']] = record['pmid']
            else:
                pmid_dict[record['pmcid']] = None
        return pmid_dict
    else:
        return None

# Example usage
pmcid_list = ["PMC1234567", "PMC7654321"]
pmid_dict = pmc_to_pmid(pmcid_list)

for pmcid, pmid in pmid_dict.items():
    print(f"PMCID: {pmcid} -> PMID: {pmid}")
