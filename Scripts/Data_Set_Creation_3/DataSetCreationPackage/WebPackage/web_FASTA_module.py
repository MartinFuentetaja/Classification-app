#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import requests 
import sys
import os

#####################################################################################################################################################################################################################

def getting_FASTA_file_ClinVar(gene_name, response, Data_path):
    query_key  = response.json()["esearchresult"]["querykey"]
    web_env    = response.json()["esearchresult"]["webenv"]
    url_1      = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&query_key={query_key}&WebEnv={web_env}&rettype=fasta&retmode=text"
    response_1 = requests.get(url_1)
    filename   = f"{gene_name}.fasta"
    filename   = os.path.join(Data_path, filename)
    with open(filename, 'wb') as f: 
        f.write(response_1.content)
    return filename


def getting_FASTA_file_UniRef(gene_name, response, Data_path):
    gene_name = gene_name.upper()
    filename   = f"{gene_name}.fasta"
    filename   = os.path.join(Data_path, filename)
    with open(filename, 'wb') as f: 
        f.write(response.content)
    return filename


def getting_FASTA_file_UniRef_2(gene_name, primaryAccession, Data_path):
    gene_name = gene_name.upper()
    FASTA_url = f"https://rest.uniprot.org/uniprotkb/{primaryAccession}.fasta"
    response  = requests.get(FASTA_url)
    filename  = f"{gene_name}.fasta"
    filename  = os.path.join(Data_path, filename)
    with open(filename, 'wb') as f:
        f.write(response.content)
    return filename
