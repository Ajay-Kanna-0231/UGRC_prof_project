from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

Entrez.email = "relaxchumma@gmail.com"

abstract = []
title = []
journal = []
year = []
month = []
language = []
doi = []
valid_pmids = []

def fetch(pm_id):
    try:
        handle = Entrez.efetch(db="pubmed", id=pm_id, rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()
    except Exception as e:
        print(f"Error fetching data for PM ID {pm_id}: {e}")
        return False
    
    root = ET.fromstring(xml_data)
    
    abstract_text = root.find('.//AbstractText')
    if abstract_text is not None and abstract_text.text is not None:
        abstract.append(abstract_text.text.strip())
    else:
        abstract.append(None)

    article_title = root.find('.//ArticleTitle')
    if article_title is not None and article_title.text is not None:
        title.append(article_title.text.strip())
    else:
        title.append(None)

    journal_title = root.find('.//Journal/Title')
    if journal_title is not None and journal_title.text is not None:
        journal.append(journal_title.text.strip())
    else:
        journal.append(None)

    language_tag = root.find('.//Language')
    if language_tag is not None and language_tag.text is not None:
        language.append(language_tag.text.strip())
    else:
        language.append(None)

    doi_found = None
    for eloc in root.findall('.//ELocationID'):
        if eloc.get('EIdType') == 'doi':
            doi_found = eloc.text.strip()
            break
    
    if doi_found is None:
        for article_id in root.findall('.//ArticleId'):
            if article_id.get('IdType') == 'doi':
                doi_found = article_id.text.strip()
                break
    
    doi.append(doi_found)

    pub_date = root.find('.//PubDate')
    if pub_date is not None and pub_date.text is not None:
        year_text = pub_date.find('Year')
        month_text = pub_date.find('Month')
        if year_text is not None and year_text.text is not None:
            year.append(year_text.text.strip())
        else:
            year.append(None)
        if month_text is not None and month_text.text is not None:
            month.append(month_text.strip())
        else:
            month.append(None)
    else:
        year.append(None)
        month.append(None)
    
    return True


'''start_pm_id = 38800000
end_pm_id = 38800100

for pm_id in range(start_pm_id, end_pm_id):
    if fetch(pm_id):
        valid_pmids.append(pm_id)'''

pm_id_list = []
for i in pm_id_list:


df = pd.DataFrame(list(zip(valid_pmids, abstract, title, journal, language, year, month, doi)),
                  columns=['pmid', 'abstract', 'title', 'journal', 'language', 'year', 'month', 'doi'])

df['full_text_urls'] = df['doi'].apply(lambda x: f"https://doi.org/{x}" if x else None)

print(df)




import json
from huggingface_hub import login
login("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")

output_file = r"C:\Users\Ajay Kanna\Desktop\UGRC\extract abstract\dump.json"
df.to_json(output_file, orient='records', lines=True)
print("Done converting to json and dumped to dump.json")

from huggingface_hub import create_repo, HfApi
api = HfApi()

repo_name = "journal-retrieval-code"
user = "data-is-everything"
try:
    # Check if the repository already exists
    api.repo_info(repo_id=f"{user}/{repo_name}")
    print(f"Repository '{repo_name}' already exists.")
    repo_url = f"https://huggingface.co/{user}/{repo_name}"
    
except RepositoryNotFoundError:
    # Create the repository if it does not exist
    repo_url = create_repo(repo_name, private=False)
    print(f"Repository created at: {repo_url}")
          
print(f"Repository URL: {repo_url}")

from huggingface_hub import upload_file
local_file_path = r"C:\Users\Ajay Kanna\Desktop\UGRC\extract abstract\dump.json"
path_in_repo = "dataset.json"  # The path where the file will be stored in the repo
repo_id = "data-is-everything/journal-retrieval-code"

upload_file(path_or_fileobj=local_file_path, path_in_repo=path_in_repo, repo_id=repo_id, token="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")
