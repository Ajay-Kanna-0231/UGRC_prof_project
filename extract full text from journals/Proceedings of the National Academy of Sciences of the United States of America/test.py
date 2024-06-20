from Bio import Entrez
from xml.etree import ElementTree as ET

def check_pmc_ids(pm_id):
    Entrez.email = "relaxchumma@gmail.com"
    
    handle = Entrez.efetch(db="pubmed", id = pm_id, rettype="full", retmode="xml")
    xml_data = handle.read()
    handle.close()
    root = ET.fromstring(xml_data)
    print(ET.tostring(root, encoding='utf8').decode('utf8'))
    pmc_dict = {}
    print("------------------------------------------------------------------------")
    pubmed_articles = root.findall(".//PubmedArticle")
    print(f"Found {len(pubmed_articles)} PubmedArticles")
    for article in root.findall(".//PubmedArticle"):
        pmid = article.find(".//PMID").text
        pmc_id = None
        for article_id in article.findall(".//ArticleId"):
            if article_id.attrib.get("IdType") == "pmc":
                pmc_id = article_id.text
                break
        pmc_dict[pmid] = pmc_id is not None
    
    return pmc_dict

print(check_pmc_ids(10377390))