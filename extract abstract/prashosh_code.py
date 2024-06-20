import numpy as np
import pandas as pd
from Bio import Entrez
import os
import json
import time
import requests
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.chrome.service import Service

Entrez.email = "ajaykannaalwar@gmail.com"  # replace with a valid email id
api_key = "794318879572f8fc0d465d5d2aa72f4f7808"  # replace with a valid api key

# service = Service()
# options = webdriver.ChromeOptions()
# driver = webdriver.Chrome(service=service, options=options)
driver = webdriver.Chrome()

def pubmed_id_list(start_id: int, end_id: int) -> list[str]:
    """
    Gets the starting PMID and end PMID.

    Args:
        start_id (int): The starting PMID.
        end_id (int): The ending PMID.

    Returns:
        list[str]: list of PMIDs where each element is a str
    """
    return list(str(i) for i in range(start_id, end_id + 1))


def post_pubmed_ids(id_list: list[str]):
    """
    Uploads a list of PMIDs to Entrez history server.

    Args:
        id_list (str[list]): The list of PMIDs need to be posted.

    Returns:
        string: WebEnv for the epost
        string: QueryKey for the epost
    """
    handle = Entrez.epost(db="pubmed", id=",".join(id_list), api_key=api_key)
    result = Entrez.read(handle)
    handle.close()
    return result["WebEnv"], result["QueryKey"]


def fetch_details(webenv, query_key, retstart, retmax=500):
    """
    Fetches the details of the papers for the given PMIDs.

    Args:
        webenv (str): The WebEnv of the EPost performed.
        query_key (str): The QueryKey of the EPost performed.
        retstart (int): the starting index for fetching.

    Returns:
        Bio.Entrez.Parser.DictionaryElement: The parsed XML results of the search in a dictionary-like object.
    """
    handle = Entrez.efetch(db="pubmed", retmode="xml", webenv=webenv, query_key=query_key, retstart=retstart, retmax=retmax, api_key=api_key)
    records = Entrez.read(handle)
    handle.close()
    return records


def get_doi(article_id_list):
    for id_element in article_id_list:
        if id_element.attributes['IdType'] == 'doi':
            return str(id_element)
    return None


def get_fulltext_link(pubmed_id, doi):
    url = f"https://pubmed.ncbi.nlm.nih.gov/{int(pubmed_id)}/"

    # Send an HTTP GET request to the website
    response = requests.get(url)

    if response.status_code == 200:
        # Parse the HTML content of the page with BeautifulSoup
        soup = BeautifulSoup(response.content, 'html.parser')

        div_element = soup.find('div', class_='full-text-links-list')
        if div_element:
            # Find all a tags inside the div element
            anchor_elements = div_element.find_all('a')

            links = []
            pmc_link = ''
            for anchor in anchor_elements:
                href = anchor.get('href')
                if "https://www.ncbi.nlm.nih.gov/pmc" in href:
                    pmc_link = href
                else:
                    links.append(href)

            if pmc_link:
                return pmc_link
            else:
                return links[0]
        else:
            if doi:
                return "https://doi.org/" + doi
            else:
                return None
    else:
        return None


def download_image(img_url, filename):
    response = requests.get(img_url)
    if response.status_code == 200:
        with open(filename, 'wb') as f:
            f.write(response.content)
    else:
        pass


def scrape_fulltext_and_images(fulltext_link, output_folder):
    print(fulltext_link)
    if not fulltext_link.startswith("http"):
        fulltext_link = "https:/" + fulltext_link

    driver.get(fulltext_link)
    html = driver.page_source
    print(html)
    soup = BeautifulSoup(html, 'html.parser')
    # article = soup.find('article')
    elements = soup.find_all(['p', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6'])
    if elements:
        full_text = []
        # Print the tag name and text content of each element in the correct order
        for element in elements:
            full_text.append(element.text.strip())
        print_text = "\n".join(full_text)
        return print_text
    else:    
        full_text = "Not Available"
    # if article:
    #     full_text = article.get_text(separator="\n", strip=True)

        # image_tags = article.find_all('img')

        # os.makedirs(output_folder, exist_ok=True)

        # for i, img in enumerate(image_tags):
        #     img_url = img.get('src')
        #     if img_url:
        #         if not img_url.startswith("http"):
        #             img_url = "https:/" + img_url
        #         image_filename = os.path.join(output_folder, f"image_{i}.jpg")
        #         download_image(img_url, image_filename)
    return full_text


def retrieve_articles(start_pubmedid, end_pubmedid, batch_size=10000, rate_limit=10):
    # pubmed_ids = pubmed_id_list(start_pubmedid, end_pubmedid)
    
    pubmed_ids = ['38836405', '38834327', '38825318', '38343800', '38822254']
    for i in range(0, len(pubmed_ids), batch_size):
        batch_ids = pubmed_ids[i:i + batch_size]
        webenv, query_key = post_pubmed_ids(batch_ids)
        retstart = 0
        while retstart < len(batch_ids):
            try:
                papers = fetch_details(webenv, query_key, retstart, retmax=500)
                for paper in papers['PubmedArticle']:
                    pmid = paper['MedlineCitation']['PMID']
                    print(pmid)
                    article_id_list = paper['PubmedData']['ArticleIdList']
                    doi = get_doi(article_id_list)
                    title = paper['MedlineCitation']['Article'].get('ArticleTitle', 'Not Available')
                    abstract_list = paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', ['Not Available'])
                    abstract = abstract_list[0] if abstract_list else 'Not Available'
                    
                    fulltext = "Not Available"  # Initialize fulltext to avoid usage before assignment
                    full_text_link = get_fulltext_link(pmid, doi)
                    output_folder = os.path.join("Text mining", "data", pmid)
                    if full_text_link:
                        fulltext = scrape_fulltext_and_images(full_text_link, output_folder)

                    details = f"PMID :  {pmid}\n\nTitle :  {title}\n\nAbstract :  {abstract}\n\nFull_Text :  {fulltext}"
                    txt_filename = os.path.join(output_folder, f"{pmid}_text.txt")
                    os.makedirs(output_folder, exist_ok=True)
                    with open(txt_filename, 'w', encoding='utf-8') as txt_file:
                        txt_file.write(details)

                retstart += 500
                time.sleep(1 / rate_limit)  # to ensure rate limit isn't exceeded  
            except Exception as e:
                print(e)
                print(f"Error in chunk {i / batch_size} starting at index {retstart}")
                retstart += 500
                time.sleep(1 / rate_limit)
                continue


retrieve_articles(start_pubmedid=1, end_pubmedid=2, batch_size=1, rate_limit=10)
driver.close()