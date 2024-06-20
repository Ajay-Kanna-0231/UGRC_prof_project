import threading
from queue import Queue
from selenium import webdriver
from bs4 import BeautifulSoup
import time
import re

import pickle
with open(r"C:\Users\Ajay Kanna\Desktop\UGRC\extract abstract\Lancet\store_uniques.pkl", 'rb') as f:
    urls = pickle.load(f)
print(len(urls))

class ScrapeThread(threading.Thread):
    def __init__(self, queue, semaphore):
        threading.Thread.__init__(self)
        self.queue = queue
        self.semaphore = semaphore

    def run(self):
        while True:
            url, j = self.queue.get()
            if url is None:
                break
            driver = webdriver.Chrome()
            try:
                driver.get(url)
                time.sleep(4)
                
                html = driver.page_source
                soup = BeautifulSoup(html, 'lxml')

                # Find the article tag
                article = soup.find('article')

                # Find all p and header (h1 to h6) tags inside the article tag
                if article:
                    elements = article.find_all(['p', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6'])

                    # Collect the text content of each element in the correct order
                    full_text = [element.text.strip() for element in elements]
                    full_text = "\n".join(full_text)

                    full_text = filter_text(full_text)

                    # Join the collected text into a single string with newlines
                    if full_text.strip():
                        print(f"Text extracted for URL {j}.")
                        with open(f'C:\\Users\\Ajay Kanna\\Desktop\\UGRC\\extract abstract\\Lancet\\full_text_stored\\extracted_text_{j}.txt', 'w', encoding='utf-8') as file:
                            file.write(full_text)
                else:
                    print(f"No article tag found for URL {j}.")
            finally:
                driver.quit()
                self.queue.task_done()
                self.semaphore.release()

def filter_text(text):
    # Remove references like (6) or [6] or References (6)
    text = re.sub(r'\bReferences? \(\d+\)', '', text)
    text = re.sub(r'\bReferences\b', '', text)
    # Remove "Cited by" with any number
    text = re.sub(r'\bCited by \(\d+\)', '', text)
    text = re.sub(r'\bCited by \d+', '', text)
    text = re.sub(r'\bCited by\b', '', text)

    if len(text.split()) < 50:
        text = ""
    
    return text

def main(urls):
    queue = Queue()
    semaphore = threading.Semaphore(4)  # Limit the number of concurrent threads to 4

    threads = []
    for i in range(len(urls)):
        t = ScrapeThread(queue, semaphore)
        t.start()
        threads.append(t)

    for j, url in enumerate(urls, start=1):
        semaphore.acquire()
        queue.put((url, j))

    queue.join()

    for i in range(len(threads)):
        queue.put(None)
    for t in threads:
        t.join()

if __name__ == "__main__":
    main(urls)
