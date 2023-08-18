#####################################################################################################################################################################################################################
#En este script vamos a tratar de extraer las tablas de https://databases.lovd.nl/shared/variants/ para ello vamos a hacer uso de BeautifulSoup4, pandas y requests

import pandas as pd
from bs4 import BeautifulSoup
import requests
import re
import math
import os

#####################################################################################################################################################################################################################

def scraping_lovd3_info(gene_name, response, Data_path):
    headers, length = scrap_header_and_size_lovd3(gene_name, response)
    mydata = pd.DataFrame(columns = headers)
    #Now we are looking for how many variants are distributed arround pages.
    page = math.ceil(float(length) / 1000)
    for i in range(1, page + 1):
        url = f"https://databases.lovd.nl/shared/variants/{gene_name}?search_var_status=%3D%22Marked%22%7C%3D%22Public%22&page_size=1000&page={i}"
        page = requests.get(url)
        soup = BeautifulSoup(page.text, 'lxml')
        table1 = soup.find('table', id=f"viewlistTable_CustomVL_VOT_VOG_{gene_name}")
        for j in table1.find_all('tr')[1:]:
            row_data = j.find_all('td')
            row = [i.text for i in row_data]
            length = len(mydata)
            mydata.loc[length] = row
    filename = f"{gene_name}_lovd3.txt"
    filename = os.path.join(Data_path, filename)
    mydata.to_csv(filename, sep = "\t", index = False)
    return filename

def scrap_header_and_size_lovd3(gene_name, response):
    soup = BeautifulSoup(response.text, 'lxml')
    table1 = soup.find('table', id=f"viewlistTable_CustomVL_VOT_VOG_{gene_name}")
    num_entries = soup.find('span', id = f'viewlistPageSplitText_CustomVL_VOT_VOG_{gene_name}')
    length = num_entries.text.strip().split()[0]
    #Let us obtain table headers and create a DataFrame
    headers = []
    for i in table1.find_all('th'):
        title = i.text.strip()
        headers.append(title)
    headers = [el.replace('\xa0',' ') for el in headers]
    return headers, length

