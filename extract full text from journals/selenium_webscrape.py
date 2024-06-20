from selenium import webdriver
from selenium.webdriver.common.by import By

driver = webdriver.Chrome()
#driver.get("https://academic.oup.com/cid/article-lookup/doi/10.1093/cid/ciae291")
driver.get("https://bmjopen.bmj.com/content/14/6/e080132.long")
#print(driver)
#print(driver.title)

headers = driver.find_elements(By.XPATH, "//h1 | //h2 | //h3 | //h4 | //h5 | //h6")
header_texts = [header.text for header in headers]

paragraphs = driver.find_elements(By.TAG_NAME, "p")
paragraph_texts = [paragraph.text for paragraph in paragraphs]

# Print extracted headers and paragraphs
print("Headers:")
for header in header_texts:
    print(header)

print("\nParagraphs:")
for paragraph in paragraph_texts:
    print(paragraph)

# Close the driver
driver.quit()