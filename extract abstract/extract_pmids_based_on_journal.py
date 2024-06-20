from Bio import Entrez
import time

# Set your email here
Entrez.email = "relaxchumma@gmail.com"

def fetch_pmids(journal_name, retmax=1000):
    # E-utilities search query
    search_query = f'"{journal_name}"[journal]'
    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList'], int(record['Count'])

def fetch_all_pmids(journal_name):
    pmids = []
    retmax = 10000  # Number of results to fetch at a time
    retstart = 0  # Starting position of the results
    total_count = fetch_pmids(journal_name, retmax=1)[1]
    
    while retstart < total_count:
        batch_pmids, _ = fetch_pmids(journal_name, retmax=retmax)
        pmids.extend(batch_pmids)
        retstart += retmax
        print(f"Retrieved {retstart}/{total_count} PMIDs")
        time.sleep(0.5)  # To avoid hitting the NCBI server too hard
        
    return pmids

# Example usage
journal_name = "PLOS ONE"
all_pmids = fetch_all_pmids(journal_name)
print(f"Total PMIDs retrieved: {len(all_pmids)}")

# Save PMIDs to a file
with open('plos_one_pmids.txt', 'w') as f:
    for pmid in all_pmids:
        f.write(f"{pmid}\n")
