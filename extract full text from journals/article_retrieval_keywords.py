from Bio import Entrez

def search(query):
    Entrez.email = 'ajaykannaalwar@gmail.com'
    handle = Entrez.esearch(db='pubmed', sort='relevance', retmax='250000', retmode='xml', term=query)
    results = Entrez.read(handle)
    handle.close()  # It's a good practice to close the handle after reading the results.
    return results

studies = search('COVID-19')
studiesIdList = studies['IdList']

def fetch_details(id_list):
  ids = ','.join(id_list)
  Entrez.email = 'ajaykannaalwar@gmail.com'
  handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
  results = Entrez.read(handle)
  return results

import pandas as pd

# List initializations
title_list = []
abstract_list = []
journal_list = []
language_list = []
pubdate_year_list = []
pubdate_month_list = []

# Function to fetch article details; assumed to be defined elsewhere
studies = fetch_details(studiesIdList)  # This line seems misplaced and is commented out.

chunk_size = 10000

# Assuming studiesIdList is defined somewhere earlier in your code.
for chunk_i in range(0, len(studiesIdList), chunk_size):
    chunk = studiesIdList[chunk_i:chunk_i + chunk_size]
    papers = fetch_details(chunk)  # Fetching details for the current chunk

    for i, paper in enumerate(papers['PubmedArticle']):
        # Title
        title_list.append(paper['MedlineCitation']['Article']['ArticleTitle'])

        # Abstract
        try:
            abstract_list.append(paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
        except:
            abstract_list.append('No Abstract')

        # Journal Title
        journal_list.append(paper['MedlineCitation']['Article']['Journal']['Title'])

        # Language
        language_list.append(paper['MedlineCitation']['Article']['Language'][0])

        # Publication Year
        try:
            pubdate_year_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
        except:
            pubdate_year_list.append('No Data')

        # Publication Month
        try:
            pubdate_month_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month'])
        except:
            pubdate_month_list.append('No Data')

# Creating a DataFrame from the accumulated data
df = pd.DataFrame(list(zip(
    title_list, abstract_list, journal_list, language_list, pubdate_year_list, pubdate_month_list
)),
columns=[
    'Title', 'Abstract', 'Journal', 'Language', 'Year', 'Month'
])

# Display the DataFrame or save to a file as needed
print(df)
