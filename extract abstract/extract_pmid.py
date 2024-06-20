from Bio import Entrez
import time

# Set your email here
Entrez.email = "relaxchumma@gmail.com"

def fetch_pmids(journal_name, retmax=10000):
    # E-utilities search query
    search_query = '("Proceedings of the National Academy of Sciences of the United States of America"[Journal])'
    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList'], int(record['Count'])

def fetch_all_pmids(journal_name):
    pmids = []
    retmax = 10000  # Number of results to fetch at a time
    total_count = fetch_pmids(journal_name, retmax=10000)[1]
    tmp = fetch_pmids(journal_name, retmax=10000)[0]
    print(total_count)
    print(len(tmp))

fetch_all_pmids("PLOS ONE")