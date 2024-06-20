import pprint
from metapub import FindIt
import os
import time

works = 0
unavailable = 0
wrong = 0


def generate_pmids(end_index):
    return list(str(i) for i in range(1, end_index + 1))


def process_pmids(pmids):
    global works, unavailable, wrong
    count = 0
    for pmid in pmids:
        if count%10 == 0:
            print(count)
        count +=1
        print(pmid)
        time.sleep(0.4)
        try:
            src = FindIt(pmid)
            if src.url == None:
                unavailable +=1
                with open('exception.txt', 'a') as exception_file:
                    exception_file.write(f"PMID: {pmid}\n")
            else:
                works += 1
                with open('works.txt', 'a') as works_file:
                    works_file.write(f"PMID: {pmid}, DOI: {src.doi}, URL: {src.url}\n")
        except:
            wrong += 1
            with open('wrongid.txt', 'a') as exception_file:
                    exception_file.write(f"PMID: {pmid}\n")

end_index = 1000
pmids = generate_pmids(end_index)
process_pmids(pmids)
print("total ids given : ",len(pmids))
print("pdf yes : ",works)
print("pdf no : ",unavailable)
print("wrong_ids : ",wrong)