#############################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd
import numpy as np
import os, shutil
import datetime
import hashlib
import requests

#####################################################################################################################################################################################################################

def obtener_hash_archivo(nombre_archivo):
    # Crea un objeto de hash SHA-256
    hash_sha256 = hashlib.sha256()

    # Abre el archivo en modo binario y lee los datos en bloques
    with open(nombre_archivo, 'rb') as archivo:
        # Lee el primer bloque de datos
        bloque = archivo.read(4096)

        # Itera mientras haya datos en el bloque
        while bloque:
            # Actualiza el hash con el bloque de datos
            hash_sha256.update(bloque)

            # Lee el siguiente bloque de datos
            bloque = archivo.read(4096)

    # Devuelve el hash en formato hexadecimal
    return hash_sha256.hexdigest()

def controlling_released_database(old_variant_summary_file):
    pwd      = os.getcwd()
    old_hash = obtener_hash_archivo(old_variant_summary_file)
    url      = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    response = requests.get(url, stream=True)
    print(f"ClinVar DataBase status: {response.status_code}")
    new_variant_summary_file = "variant_summary(1).txt.gz"
    path = os.path.join(pwd, new_variant_summary_file)
    with open(new_variant_summary_file, 'wb') as f:
        f.write(response.content)
    new_hash = obtener_hash_archivo(new_variant_summary_file)
    if(old_hash == new_hash):
        os.remove(path)
        print("There are not modifications!")
    else:
        path_new         = os.path.join(pwd, "variant_summary.txt.gz")
        old_variant_path = os.path.dirname(old_variant_summary_file)
        os.remove(old_variant_summary_file)
        os.rename(path, path_new)
        shutil.move(path_new, old_variant_path)
        print("DOWNLOADED!")

"""
Estas función se encarga de extraer la información relativa de un gen en específico de la base de datos: "variant_summary.txt.gz" descargada de ncbi. De esta manera,
la función extrae todos los datos del archivo y acto seguido busca toda la información relevante al input introducido, es decir, al nombre del gen. Por otra parte,
hemos reducido el número de columnas a usecols = ["#AlleleID", "Name", "GeneSymbol", "ClinicalSignificance", "PhenotypeList", "ReviewStatus"] dado que son las que vamos a utilizar,
y por lo tanto, por un ejercicio de eficacia hemos evitado las demás. Una vez tenemos la información relativa al gen específicado, la función modificará el nombre de las columnas
y añadirá una nueva ("Protein change") con el fin de que este DataFrame tenga la misma estrucutura que el obtenido introduciendo nosotros a mano el archivo de clinvar, o bien el obtenido
de internet.
"""
def ClinVar_file_from_database(gene_name, Data_path):
    file = os.path.join(Data_path, "variant_summary.txt.gz")
    df   = pd.read_csv(file, sep = "\t", compression = "gzip", header = 0, dtype = {"Chromosome" : str}, usecols = ["#AlleleID", "Name", "GeneSymbol", "ClinicalSignificance", "PhenotypeList", "ReviewStatus"])
    df_gene   = df[df["GeneSymbol"].str.upper() == gene_name.upper()].drop_duplicates(subset = ["#AlleleID"], keep = "first", ignore_index = True)
    df_gene["Protein change"] = np.nan
    df_gene = df_gene.rename(columns={"GeneSymbol": "Gene(s)", "ClinicalSignificance": "Clinical significance (Last reviewed)", "ReviewStatus" : "Review status" ,"PhenotypeList" : "Condition(s)"})
    df_gene = df_gene[["Name", "Gene(s)", "Protein change", "Condition(s)", "Clinical significance (Last reviewed)", "Review status"]]
    file_name = f"clinvar_result_{gene_name}.txt"
    file_name = os.path.join(Data_path, file_name)
    df_gene.to_csv(file_name, sep = "\t", index = False) 
    return file_name

def date_control():
    # Obtén la fecha actual
    fecha_actual = datetime.datetime.now()

    # Obtén el día de la semana de la fecha actual (0 es lunes, 3 es jueves)
    dia_semana = fecha_actual.weekday()

    # Verifica si es jueves y si es el primer jueves del mes
    if dia_semana == 3 and fecha_actual.day <= 7:
        print("")
        print("It is the fist thursday of the month. Thus, ClinVar database is released.")
        old_variant_summary_file = "/home/einstein/martin/git_enviroment/Classification_app/Data/variant_summary.txt.gz"
        controlling_released_database(old_variant_summary_file)
