#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import requests 
import pandas as pd
import os

#####################################################################################################################################################################################################################


def obtaining_ClinVar_file(gene_name, count, Data_path):
    starting = 0
    variation_name = []
    gene_list      = []
    protein_change = []
    description    = []
    condition_tot  = []
    review_status  = []
    while (int(count) >= starting):
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}[gene]+AND+single_gene[prop]&sort=position&retmax=500&retstart={starting}&retmode=json"
        response = requests.get(url)
        gene_id = response.json()["esearchresult"]["idlist"]
        gen = ",".join(gene_id)
        starting += 500
        url_1 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={gen}&rettype=uilist&retmode=json"
        response_1 = requests.get(url_1)
        for i in range(len(gene_id)):
            condition      = []
            variation_name.append(response_1.json()["result"][gene_id[i]]["title"])
            gene_list.append(response_1.json()["result"][gene_id[i]]["genes"][0]["symbol"])
            protein_change.append(response_1.json()["result"][gene_id[i]]["protein_change"])
            description.append(response_1.json()["result"][gene_id[i]]["clinical_significance"]["description"] + f"({response_1.json()['result'][gene_id[i]]['clinical_significance']['last_evaluated']})")
            for j in range(len(response_1.json()["result"][gene_id[i]]["trait_set"])):
                condition.append(response_1.json()["result"][gene_id[i]]["trait_set"][j]["trait_name"])
            condition_tot.append("|".join(condition))
            review_status.append(response_1.json()["result"][gene_id[i]]["clinical_significance"]["review_status"])
    dictionary = {'Name' : variation_name, 'Gene(s)' : gene_list, 'Protein change' : protein_change, "Condition(s)" : condition_tot, "Clinical significance (Last reviewed)" : description, "Review status" : review_status}
    df = pd.DataFrame(dictionary)
    file_name = f"clinvar_result_{gene_name}.txt"
    file_name = os.path.join(Data_path, file_name)
    df.to_csv(file_name, sep = "\t", index = False) 
    return file_name

