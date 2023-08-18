import requests
import pandas as pd
import sys


def extract_id_ensembl(gene_name, length):
    gene_name = gene_name.upper()
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1;content-type=application/json"
    response = requests.get(url)
    if (response.status_code == 200 and len(response.json()) > 1):
        transcript = extracting_info(response, length)
    elif (response.status_code != 200):
        print(f"Ensemble response status: {response.status_code}")
        sys.exit(f"There has been a problem extracting information from ensembl.org. Please check that the gene_name is correct. {response.status_code}")
    elif (len(response.json()) == 1):
        print(f"Ensemble information status: {response.json()}")
        sys.exit(f"There has been a problem extracting information from ensembl.org. Please check that the gene_name is correct. {response.json()}")
    return transcript

def extracting_info(response, length):
    transcript = response.json()["Transcript"]
    for i in range(len(transcript)):
        if (transcript[i]["biotype"] == "protein_coding"):
            try:
                returned_length = int(transcript[i]["Translation"]["length"])
                if (int(length) == returned_length):
                    return transcript[i]["id"] + "." + str(transcript[i]["version"])
            except KeyError:
                sys.exit(f"There does not exist any gene with the indicated length: {length}")
    sys.exit(f"The length for the indicated gene corresponds to a not protein coding gene.")
    
