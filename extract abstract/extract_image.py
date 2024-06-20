from selenium import webdriver
from bs4 import BeautifulSoup
import requests
import os

def get_all_images(full_url, download=False):
    # Initialize the Chrome driver
    driver = webdriver.Chrome()
    driver.get(full_url)

    # Get the page source and parse it with BeautifulSoup
    html = driver.page_source
    soup = BeautifulSoup(html, 'html.parser')

    # Find all image elements
    img_tags = soup.find_all('img')

    # Extract the URLs of the images, excluding .svg files
    img_urls = []
    for img in img_tags:
        if 'src' in img.attrs:
            src = img['src']
            if not src.endswith('.svg'):
                if src.startswith('http'):
                    img_urls.append(src)
                else:
                    # Handle relative URL by making it absolute
                    img_urls.append(requests.compat.urljoin(full_url, src))

    # Print the URLs
    for img_url in img_urls:
        print(img_url)

    # Optionally, download the images
    if download:
        if not os.path.exists('images'):
            os.makedirs('images')

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }

        for i, img_url in enumerate(img_urls):
            try:
                response = requests.get(img_url, headers=headers)
                response.raise_for_status()

                img_data = response.content
                with open(f'images/image_{i+1}.jpg', 'wb') as img_file:
                    img_file.write(img_data)
            except requests.HTTPError as e:
                print(f"HTTP error occurred when downloading {img_url}: {e}")
            except requests.RequestException as e:
                print(f"Request error occurred for {img_url}: {e}")
            except Exception as e:
                print(f"An error occurred when handling {img_url}: {e}")
    
    # Close the driver
    driver.quit()

# URL to process
full_url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11127088/"  # Replace with your target URL

# Call the function
get_all_images(full_url, download=True)  # Set download to True to download images
