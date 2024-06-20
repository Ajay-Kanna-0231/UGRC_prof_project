from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
import time
import pandas as pd

driver = webdriver.Chrome()

driver.get("https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/38800010/")

html = driver.page_source
soup = BeautifulSoup(html, 'html.parser')

p_tags = soup.find_all(['p', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6'])
text = ' '.join([p.get_text() for p in p_tags])

tables = soup.find_all('table')
for table in tables:
    headers = [th.text.strip() for th in table.find_all('th')]
    rows = []
    for tr in table.find_all('tr'):
        cells = tr.find_all('td')
        if not cells:
            continue
        row = [cell.text.strip() for cell in cells]
        rows.append(row)

    df = pd.DataFrame(rows, columns=headers)
    print(df)

with open('extracted_text.txt', 'w', encoding='utf-8') as file:
    file.write(text)

driver.quit()
