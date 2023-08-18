#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""


import DataSetCreationPackage.WebPackage.web_ClinVar_module   as web_ClinVar_mod
import DataSetCreationPackage.WebPackage.web_AlphaFold_module as web_AlphaFold_mod
import DataSetCreationPackage.WebPackage.web_FASTA_module     as web_FASTA_mod
import DataSetCreationPackage.WebPackage.web_Lovd3_module     as web_Lovd3_mod
import DataSetCreationPackage.WebPackage.web_GenomAD_module   as web_GenomAD_mod
import DataSetCreationPackage.WebPackage.web_Ensembl_module   as web_ensembl_mod
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

def initializing_AlphaFold_file(response, length, gene_name, Data_path):
    primaryAccession, returned_length = web_AlphaFold_mod.obtain_primaryAccession_uniProt_modified(response, length, gene_name)
    AlphaFold_file = web_AlphaFold_mod.downloading_PDB_file_AlphaFold(primaryAccession, Data_path)
    return AlphaFold_file, primaryAccession, returned_length

def controlling_status(url):
    response    = requests.get(url)
    status_code = response.status_code
    return response, status_code

def generating_url(gene_name, length):
    gene_name      = gene_name.upper()
    length_initial = int(length) - 200
    length_final   = int(length) + 200
    #Extraemos el transcriptID para hacer la busqueda en GenomAD. Primero debemos hacer un split en "." y coger el primero elemento.
    transcript_id  = web_ensembl_mod.extract_id_ensembl(gene_name, int(length)).split(".")[0]
    UniProt_url    = f"https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28gene%3A{gene_name}%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28length%3A%5B{length_initial}%20TO%20{length_final}%5D%29&size=500"
    ClinVar_url    = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}[gene]+AND+single_gene[prop]&retmax=20&retmode=json"
    Lovd3_url      = f"https://databases.lovd.nl/shared/variants/{gene_name}?search_var_status=%3D%22Marked%22%7C%3D%22Public%22&page_size=10&page=1"
    genomAD_V2_url           = f"https://gnomad.broadinstitute.org/transcript/{transcript_id}?dataset=gnomad_r2_1"
    genomAD_V2_control_url   = f"https://gnomad.broadinstitute.org/transcript/{transcript_id}?dataset=gnomad_r2_1_controls"
    genomAD_V2_non_neuro_url = f"https://gnomad.broadinstitute.org/transcript/{transcript_id}?dataset=gnomad_r2_1_non_neuro"
    genomAD_V3_non_V2_url    = f"https://gnomad.broadinstitute.org/transcript/{transcript_id}?dataset=gnomad_r3_non_v2"
    return UniProt_url, ClinVar_url, Lovd3_url, genomAD_V2_url, genomAD_V2_control_url, genomAD_V2_non_neuro_url, genomAD_V3_non_V2_url

def obtaining_files_web(gene_name, length, Data_Path):
    #Inicializamos cada url para la descarga de ficheros
    UniProt_url, ClinVar_url, Lovd3_url, genomAD_V2_url, genomAD_V2_control_url, genomAD_V2_non_neuro_url, genomAD_V3_non_V2_url  = generating_url(gene_name, length)
    #Obtenemos las respuestas recibidas por cada API
    response_UniProt, status_code_UniProt = controlling_status(UniProt_url)
    response_ClinVar, status_code_ClinVar = controlling_status(ClinVar_url)
    response_Lovd3, status_code_Lovd3     = controlling_status(Lovd3_url)
    response_genomAD_V2, status_code_genomAD_V2                 = controlling_status(genomAD_V2_url)
    response_genomAD_V2_control, status_code_genomAD_V2_control = controlling_status(genomAD_V2_control_url)
    response_genomAD_V2_neuro, status_code_genomAD_V2_neuro     = controlling_status(genomAD_V2_non_neuro_url)
    response_genomAD_V3, status_code_genomAD_V3                 = controlling_status(genomAD_V3_non_V2_url)
    #Esta parte se encarga de analizar únicamente una parte de la página de GenomAD. Controla si en el contenido de la página web no aparece nada.
    response_genomAD_V2         = web_GenomAD_mod.scraping_genomAD_info(genomAD_V2_url)
    response_genomAD_V2_control = web_GenomAD_mod.scraping_genomAD_info(genomAD_V2_control_url)
    response_genomAD_V2_neuro   = web_GenomAD_mod.scraping_genomAD_info(genomAD_V2_non_neuro_url)
    response_genomAD_V3         = web_GenomAD_mod.scraping_genomAD_info(genomAD_V3_non_V2_url)
    if ((status_code_UniProt == 200 and len(response_UniProt.json()["results"]) != 0) and (status_code_ClinVar == 200 and int(response_ClinVar.json()["esearchresult"]["count"]) != 0) and (status_code_Lovd3 == 200 and response_Lovd3.text != "") 
        and (status_code_genomAD_V2 == 200 and response_genomAD_V2 != "Gene not found") and (status_code_genomAD_V2_control == 200 and response_genomAD_V2_control != "Gene not found") and (status_code_genomAD_V2_neuro == 200 and response_genomAD_V2_neuro != "Gene not found") and (status_code_genomAD_V3 == 200 and response_genomAD_V3 != "Gene not found")):
        #Una vez comprobada que las respuestas son correctas, se descargan los archivos correspondientes.
        AlphaFold_file, primaryAccession, returned_length = initializing_AlphaFold_file(response_UniProt, length, gene_name, Data_Path)
        FASTA_file   = web_FASTA_mod.getting_FASTA_file_UniRef_2(gene_name, primaryAccession, Data_Path)
        ClinVar_file = web_ClinVar_mod.obtaining_ClinVar_file(gene_name, int(response_ClinVar.json()["esearchresult"]["count"]), Data_Path)
        Lovd3_file   = web_Lovd3_mod.scraping_lovd3_info(gene_name, response_Lovd3, Data_Path)
        GenomAD_V2_file           = web_GenomAD_mod.download_genomAD(genomAD_V2_url, Data_Path)
        GenomAD_V2_control_file   = web_GenomAD_mod.download_genomAD(genomAD_V2_control_url, Data_Path)
        GenomAD_V2_non_neuro_file = web_GenomAD_mod.download_genomAD(genomAD_V2_non_neuro_url, Data_Path)
        GenomAD_V3_file           = web_GenomAD_mod.download_genomAD(genomAD_V3_non_V2_url, Data_Path)
        print("DOWNLOADED!")
    elif ((status_code_UniProt != 200) or (status_code_ClinVar != 200) or (status_code_Lovd3 != 200) or (status_code_genomAD_V2 != 200) or (status_code_genomAD_V2_control != 200) or (status_code_genomAD_V2_neuro != 200) or (status_code_genomAD_V3 != 200)):
        #En caso de haber un error de comunicación con el servidor (por no estar disponible o mantenimientos) salatará el siguiente error y la ejecución finalizará
        print("")
        print(f"AlphaFold status:            {status_code_UniProt}")
        print(f"ClinVar status:              {status_code_ClinVar}")
        print(f"ClinVar status:              {status_code_Lovd3}")
        print(f"GenomAD_V2 status:           {status_code_genomAD_V2}")
        print(f"GenomAD_V2_control status:   {status_code_genomAD_V2_control}")
        print(f"GenomAD_V2_non_neuro status: {status_code_genomAD_V2_neuro}")
        print(f"GenomAD_V3 status:           {status_code_genomAD_V3}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    elif ((len(response_UniProt.json()["results"]) == 0) or (int(response_ClinVar.json()["esearchresult"]["count"]) == 0) or (response_Lovd3.text == "") or (response_genomAD_V2 == "Gene not found") or (response_genomAD_V3 == "Gene not found")):
        print("")
        print(f"AlphaFold information length:     {len(response_UniProt.json()['results'])}")
        print(f"ClinVar information length:       {int(response_ClinVar.json()['esearchresult']['count'])}")
        print(f"Lovd3 information length:         {response_Lovd3.text}")
        print(f"GenomAD_V2 information:           {response_genomAD_V2}")
        print(f"GenomAD_V2_control information:   {response_genomAD_V2_control}")
        print(f"GenomAD_V2_non_neuro information: {response_genomAD_V2_neuro}")
        print(f"GenomAD_V3 information:           {response_genomAD_V3}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    return AlphaFold_file, returned_length, FASTA_file, ClinVar_file, Lovd3_file, GenomAD_V2_file, GenomAD_V2_control_file, GenomAD_V2_non_neuro_file, GenomAD_V3_file


def obtaining_files_web_DataBase(gene_name, length, Data_path):
    #Inicializamos cada url para la descarga de ficheros
    UniProt_url, ClinVar_url, Lovd3_url, genomAD_V2_url, genomAD_V2_control_url, genomAD_V2_non_neuro_url, genomAD_V3_non_V2_url = generating_url(gene_name, length)
    #Obtenemos las respuestas recibidas por cada API
    response_UniProt, status_code_UniProt = controlling_status(UniProt_url)
    response_Lovd3, status_code_Lovd3     = controlling_status(Lovd3_url)
    response_genomAD_V2, status_code_genomAD_V2 = controlling_status(genomAD_V2_url)
    response_genomAD_V2_control, status_code_genomAD_V2_control = controlling_status(genomAD_V2_control_url)
    response_genomAD_V2_neuro, status_code_genomAD_V2_neuro     = controlling_status(genomAD_V2_non_neuro_url)
    response_genomAD_V3, status_code_genomAD_V3 = controlling_status(genomAD_V3_non_V2_url)
    #Esta parte se encarga de analizar únicamente una parte de la página de GenomAD. Controla si en el contenido de la página web no aparece nada.
    response_genomAD_V2         = web_GenomAD_mod.scraping_genomAD_info(genomAD_V2_url)
    response_genomAD_V2_control = web_GenomAD_mod.scraping_genomAD_info(genomAD_V2_control_url)
    response_genomAD_V2_neuro   = web_GenomAD_mod.scraping_genomAD_info(genomAD_V2_non_neuro_url)
    response_genomAD_V3         = web_GenomAD_mod.scraping_genomAD_info(genomAD_V3_non_V2_url)
    if ((status_code_UniProt == 200 and len(response_UniProt.json()["results"]) != 0) and (status_code_Lovd3 == 200 and response_Lovd3.text != "") and (status_code_genomAD_V2 == 200 and response_genomAD_V2 != "Gene not found") and (status_code_genomAD_V2_control == 200 and response_genomAD_V2_control != "Gene not found") 
        and (status_code_genomAD_V2_neuro == 200 and response_genomAD_V2_neuro != "Gene not found") and (status_code_genomAD_V3 == 200 and response_genomAD_V3 != "Gene not found")):
        #Una vez comprobada que las respuestas son correctas, se descargan los archivos correspondientes.
        AlphaFold_file, primaryAccession, returned_length = initializing_AlphaFold_file(response_UniProt, length, gene_name, Data_path)
        FASTA_file   = web_FASTA_mod.getting_FASTA_file_UniRef_2(gene_name, primaryAccession, Data_path)
        ClinVar_file = ClinVar_data_mod.ClinVar_file_from_database(gene_name, Data_path)
        Lovd3_file   = web_Lovd3_mod.scraping_lovd3_info(gene_name, response_Lovd3, Data_path)
        GenomAD_V2_file           = web_GenomAD_mod.download_genomAD(genomAD_V2_url, Data_path)
        GenomAD_V2_control_file   = web_GenomAD_mod.download_genomAD(genomAD_V2_control_url, Data_path)
        GenomAD_V2_non_neuro_file = web_GenomAD_mod.download_genomAD(genomAD_V2_non_neuro_url, Data_path)
        GenomAD_V3_file           = web_GenomAD_mod.download_genomAD(genomAD_V3_non_V2_url, Data_path)
        print("DOWNLOADED!")
    elif ((status_code_UniProt != 200) or (status_code_Lovd3 != 200) or (status_code_genomAD_V2 != 200) or (status_code_genomAD_V3 != 200)):
        #En caso de haber un error de comunicación con el servidor (por no estar disponible o mantenimientos) salatará el siguiente error y la ejecución finalizará
        print("")
        print(f"AlphaFold status:            {status_code_UniProt}")
        print(f"Lovd3 status:                {status_code_Lovd3}")
        print(f"GenomAD_V2 status:           {status_code_genomAD_V2}")
        print(f"GenomAD_V2_control status:   {status_code_genomAD_V2_control}")
        print(f"GenomAD_V2_non_neuro status: {status_code_genomAD_V2_neuro}")
        print(f"GenomAD_V3 status:           {status_code_genomAD_V3}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    elif ((len(response_UniProt.json()["results"]) == 0) or (response_Lovd3.text == "") or (response_genomAD_V2 == "Gene not found") or (response_genomAD_V3 == "Gene not found")):
        print("")
        print(f"AlphaFold information length: {len(response_UniProt.json()['results'])}")
        print(f"Lovd3 information length:     {response_Lovd3.text}")
        print(f"GenomAD_V2 information:       {response_genomAD_V2}")
        print(f"GenomAD_V2_control information:   {response_genomAD_V2_control}")
        print(f"GenomAD_V2_non_neuro information: {response_genomAD_V2_neuro}")
        print(f"GenomAD_V3 information:       {response_genomAD_V3}")
        sys.exit(f"There has been a problem with downloading files. Maybe, you have written gen_name badly.")
    return AlphaFold_file, returned_length, FASTA_file, ClinVar_file, Lovd3_file, GenomAD_V2_file, GenomAD_V2_control_file, GenomAD_V2_non_neuro_file, GenomAD_V3_file


