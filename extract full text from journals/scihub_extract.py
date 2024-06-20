import os
import requests
from bs4 import BeautifulSoup
import wget
from time import sleep
import pprint
from metapub import FindIt

output_folder = r'C:\Users\Ajay Kanna\Desktop\UGRC\extract abstract\Downloaded_From_scihub'

# Ensure the output folder exists
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def download_paper(doi, output_folder, pmid):
    try:
        base_url = 'https://sci-hub.se/'
        response = requests.get(base_url + doi.strip())
        soup = BeautifulSoup(response.content, 'lxml')
        content = soup.find('embed').get('src').replace('#navpanes=0&view=FitH', '').replace('//', '/')
        
        if content.startswith('/downloads'):
            pdf = 'https://sci-hub.se' + content
        elif content.startswith('/tree'):
            pdf = 'https://sci-hub.se' + content
        elif content.startswith('/uptodate'):
            pdf = 'https://sci-hub.se' + content
        else:
            pdf = 'https:/' + content

        # Save the PDF to the specified directory
        file_name = os.path.join(output_folder, doi.replace('/', '_') + '.pdf')
        wget.download(pdf, out=file_name)

        # Log the successful download
        with open('PDFs_Found.txt', 'a') as pdfs:
            pdfs.write(doi.strip() + '\t' + ',' + pmid + '\t' + pdf + '\n')

    except Exception as e:
        # Log the failed download attempt
        with open('PDFs_Not_Found.txt', 'a') as nopdfs:
            nopdfs.write(doi.strip() + '\t' + ',' +  pmid + '\n')
        print(f"Failed to download DOI {doi.strip()}: {e}")

    sleep(3)

pmids = list(str(i) for i in range(38000433, 38010000))
for pmid in pmids:
    src = FindIt(pmid)
    tmp = src.doi
    if tmp:
        download_paper(tmp, output_folder, pmid)
        print("Done", pmid)
    else:
        print("DOI doesn't exist", pmid)
