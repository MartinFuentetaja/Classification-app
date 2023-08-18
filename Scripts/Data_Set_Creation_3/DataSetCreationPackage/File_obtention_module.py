#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""


import DataSetCreationPackage.WebPackage.web_ClinVar_module   as web_ClinVar_mod
import DataSetCreationPackage.WebPackage.web_AlphaFold_module as web_AlphaFold_mod
import DataSetCreationPackage.WebPackage.web_FASTA_module     as web_FASTA_mod
import DataSetCreationPackage.clinvar_database_module         as ClinVar_data_mod
import requests
import sys

#####################################################################################################################################################################################################################

"""
Esta función recoge la respuesta producida por la función "controlling_status(gene_name, length)" y el numero de variantes a analizar ("length").
De esta forma, la función "web_mod.obtain_primaryAccession_uniProt_modified(response, length)" devolverá la manera en que UniProt clasifica el gen,
es decir, nos devolverá el número que indica la clasificación del gen, y este es necesario para descargar el archivo .pdb de AlphaFold, dado que
la estructuración de nombres de los archivos son "AF-O43525-F1-model_v4.pdb", y por lo tanto, el número a obtener es "O43525". Por otro lado, también
devuelve la longitud de variantes que ha encontrado dado que puede que no correspondan al indicado.
"""

def initializing_AlphaFold_file(response, length):
    primaryAccession, returned_length = web_AlphaFold_mod.obtain_primaryAccession_uniProt_modified(response, length)
    AlphaFold_file = web_AlphaFold_mod.downloading_PDB_file_AlphaFold(primaryAccession)
    return AlphaFold_file, returned_length

def controlling_status(url):
    response    = requests.get(url)
    status_code = response.status_code
    return response, status_code

def generating_url(gene_name, length):
    gene_name      = gene_name.upper()
    length_initial = int(length) - 200
    length_final   = int(length) + 200
    AlphaFold_url     = f"https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28gene%3A{gene_name}%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28length%3A%5B{length_initial}%20TO%20{length_final}%5D%29&size=500"
    ClinVar_url       = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}[gene]+AND+single_gene[prop]&retmax=20&retmode=json"
    FASTA_url_ClinVar = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={gene_name}[gene]&usehistory=y&retmode=json"
    FASTA_url_UniRef  = f"https://rest.uniprot.org/uniref/stream?format=fasta&query=%28{gene_name}%29%20AND%20%28identity%3A1.0%29&sort=length%20asc"
    return AlphaFold_url, ClinVar_url, FASTA_url_ClinVar, FASTA_url_UniRef

def obtaining_files_web(gene_name, length):
    AlphaFold_url, ClinVar_url, FASTA_url_ClinVar, FASTA_url_UniRef = generating_url(gene_name, length)
    response_Alpha, status_code_Alpha     = controlling_status(AlphaFold_url)
    response_ClinVar, status_code_ClinVar = controlling_status(ClinVar_url)
    response_FASTA, status_code_FASTA     = controlling_status(FASTA_url_UniRef)    
    if ((status_code_Alpha == 200 and len(response_Alpha.json()["results"]) != 0) and (status_code_ClinVar == 200 and int(response_ClinVar.json()["esearchresult"]["count"]) != 0) and (status_code_FASTA == 200 and response_FASTA.content != "")): #int(response_FASTA.json()["esearchresult"]["count"]) != 0)): #and response_FASTA.content != ""
        AlphaFold_file, returned_length = initializing_AlphaFold_file(response_Alpha, length)
        #FASTA_file                      = web_FASTA_mod.getting_FASTA_file_ClinVar(gene_name, response_FASTA)
        FASTA_file                      = web_FASTA_mod.getting_FASTA_file_UniRef(gene_name, response_FASTA)
        ClinVar_file                    = web_ClinVar_mod.obtaining_ClinVar_file(gene_name, int(response_ClinVar.json()["esearchresult"]["count"]))
    elif ((status_code_Alpha != 200) or (status_code_ClinVar != 200) or (status_code_FASTA != 200)):
        print("")
        print(f"AlphaFold status: {status_code_Alpha}")
        print(f"ClinVar status:   {status_code_ClinVar}")
        print(f"FASTA status:     {status_code_FASTA}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    elif ((len(response_Alpha.json()["results"]) == 0) or (int(response_ClinVar.json()["esearchresult"]["count"]) == 0) or  response_FASTA.content == ""):#(int(response_FASTA.json()["esearchresult"]["count"]) == 0)):
        print("")
        print(f"AlphaFold information length: {len(response_Alpha.json()['results'])}")
        print(f"ClinVar information length:   {int(response_ClinVar.json()['esearchresult']['count'])}")
        print(f"FASTA information length:     {response_FASTA.content}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    return AlphaFold_file, returned_length, FASTA_file, ClinVar_file


def obtaining_files_web_DataBase_ClinVar(gene_name, length):
    AlphaFold_url, ClinVar_url, FASTA_url_ClinVar, FASTA_url_UniRef = generating_url(gene_name, length)
    response_Alpha, status_code_Alpha     = controlling_status(AlphaFold_url)
    response_ClinVar, status_code_ClinVar = controlling_status(ClinVar_url)
    response_FASTA, status_code_FASTA     = controlling_status(FASTA_url_ClinVar)
    if ((status_code_Alpha == 200 and len(response_Alpha.json()["results"]) != 0) and (status_code_FASTA == 200 and int(response_FASTA.json()["esearchresult"]["count"]) != 0)):
        print("")
        print(f"AlphaFold status: {status_code_Alpha}")
        print(f"FASTA status:     {status_code_FASTA}")
        AlphaFold_file, returned_length = initializing_AlphaFold_file(response_Alpha, length)
        FASTA_file                      = web_FASTA_mod.getting_FASTA_file_ClinVar(gene_name, response_FASTA)
        ClinVar_file                    = ClinVar_data_mod.ClinVar_file_from_database(gene_name)
        print("DOWNLOADED!")
    elif ((status_code_Alpha != 200) or (status_code_FASTA != 200)):
        print("")
        print(f"AlphaFold status: {status_code_Alpha}")
        print(f"FASTA status:     {status_code_FASTA}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    elif ((len(response_Alpha.json()["results"]) == 0) or (int(response_FASTA.json()["esearchresult"]["count"]) == 0)):
        print("")
        print(f"AlphaFold information length: {len(response_Alpha.json()['results'])}")
        print(f"FASTA information length:     {int(response_FASTA.json()['esearchresult']['count'])}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    return AlphaFold_file, returned_length, FASTA_file, ClinVar_file

def obtaining_files_web_DataBase_UniRef(gene_name, length):
    AlphaFold_url, ClinVar_url, FASTA_url_ClinVar, FASTA_url_UniRef = generating_url(gene_name, length)
    response_Alpha, status_code_Alpha     = controlling_status(AlphaFold_url)
    response_ClinVar, status_code_ClinVar = controlling_status(ClinVar_url)
    response_FASTA, status_code_FASTA     = controlling_status(FASTA_url_UniRef)
    if ((status_code_Alpha == 200 and len(response_Alpha.json()["results"]) != 0) and (status_code_FASTA == 200 and response_FASTA.content != "")): 
        print("")
        print(f"AlphaFold status: {status_code_Alpha}")
        print(f"FASTA status:     {status_code_FASTA}")
        AlphaFold_file, returned_length = initializing_AlphaFold_file(response_Alpha, length)
        FASTA_file                      = web_FASTA_mod.getting_FASTA_file_UniRef(gene_name, response_FASTA)
        ClinVar_file                    = ClinVar_data_mod.ClinVar_file_from_database(gene_name)
        print("DOWNLOADED!")
    elif ((status_code_Alpha != 200) or (status_code_FASTA != 200)):
        print("")
        print(f"AlphaFold status: {status_code_Alpha}")
        print(f"FASTA status:     {status_code_FASTA}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    elif ((len(response_Alpha.json()["results"]) == 0) or response_FASTA.content != ""): 
        print("")
        print(f"AlphaFold information length: {len(response_Alpha.json()['results'])}")
        print(f"FASTA information length:     {int(response_FASTA.json()['esearchresult']['count'])}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    return AlphaFold_file, returned_length, FASTA_file, ClinVar_file
