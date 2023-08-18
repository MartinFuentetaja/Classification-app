#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import requests 
import os
import sys

#####################################################################################################################################################################################################################

"""
Esta es una primera versión para obtener el número de clasificación (primaryAccession) que utilizan en UniProt para un gen. De esta manera,
recoge la respuesta del API y el nombre del gen, y analiza la información obtenida. Para ello, estudiamos la lista "response.json()["results"]",
la cual, es la respuesta devuelta por el API con la información que nos interesa. Esta lista está compuesta de sublistas, que corresponden a cada gen
que ha encontrado con las caracteristicas indicadas. A su vez, estas sublistas están compuestas de diccionarios con la información para cada gen. 
De esta manera, analizamos cada respuesta que ha encontrado individualmente con un for, y luego evaluamos la componente "uniProtkbId" de cada diccionario
para cada sublista. En esta componente está indicado el nombre los diferentes resultados encontrados, pero solamente analizamos el que coincida con
el nombre del gen introducido. Una vez, haya encontrado dicho gen, obtiene el número con el cual se clasifica (primaryAccession) y también el 
número de variantes que tiene (returned_length).
"""

def obtain_primaryAccession_uniProt(response, gene_name):
    size = len(response.json()["results"])
    for i in range(size):
        if (response.json()["results"][i]["uniProtkbId"].split("_")[0] == gene_name):
            primaryAccession = response.json()["results"][i]["primaryAccession"]
            returned_length  =  response.json()["results"][i]["sequence"]["length"]
    return primaryAccession, returned_length

"""
Esta es la versión modificada de la anterior función. Lo que ocurre es que en el anterior caso solo buscamos la coincidencia entre el nombre del gen
y el nombre que le de UniProt. Por lo tanto, si nuestro gen es "KCNQ3", a pesar de que UniProt ofrezca diferentes entradas para dicho gen, la función
solo devolverá la información para el caso "KCNQ3_HUMAN". Lo que ocurre, es que podemos estar interesados en analizar un gen, pero con diferentes número
de variantes. Por lo tanto, debemos estudiar la longitud indicada. Así pues, esta función encuentra el gen con la longitud de variantes indicada, o bien
el resultado que más se acerce a dicha longitud.
"""

def obtain_primaryAccession_uniProt_modified(response, length, gene_name):
    size            = len(response.json()["results"])
    d_length_memory = 10000
    index_memory    = 0
    for i in range(size):
        returned_length  = int(response.json()["results"][i]["sequence"]["length"])
        d_length         = abs(int(length) - int(returned_length))
        if (returned_length == int(length) and response.json()["results"][i]["uniProtkbId"].split("_")[0].upper() == gene_name.upper()):
            primaryAccession = response.json()["results"][i]["primaryAccession"]
            returned_length  =  length
            return primaryAccession, returned_length
        elif (d_length <= d_length_memory):
            d_length_memory = d_length
            index_memory    = i
    primaryAccession = response.json()["results"][index_memory]["primaryAccession"]
    returned_length  = response.json()["results"][index_memory]["sequence"]["length"]
    return primaryAccession, returned_length

"""
Esta función recoge el nombre del gen a estudiar y la longitud de variantes que queremos, le lanza la petición al API, y analizamos la respuesta obtenida.

La respuestas se pueden clasificar:

200 = The request was processed successfully.
400	= Bad request. There is a problem with your input.
404	= Not found. The resource you requested doesn't exist.
410	= Gone. The resource you requested was removed.
500	= Internal server error. Most likely a temporary problem, but if the problem persists please contact us.
503	= Service not available. The server is being updated, try again later.

Esto se hace para hacer un control a la hora de ejecutar el programa. De esta manera, si la respuesta es 200 el progrma se ejecutará, y sino, no.
"""

def controlling_status(gene_name, length):
    gene_name      = gene_name.upper()
    length_initial = int(length) - 200
    length_final   = int(length) + 200
    url = f"https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28gene%3A{gene_name}%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28length%3A%5B{length_initial}%20TO%20{length_final}%5D%29&size=500"
    response    = requests.get(url)
    status_code = response.status_code
    return response, status_code

"""
Una vez sabemos que la respuesta del API es correcta y que hemos obtenido el número de clasificación de UniProt para un gen con cierta longitud de 
variantes, descargamos la información de AlphaFold y la guardamos en un archivo .pdb.
"""

def downloading_PDB_file_AlphaFold(primaryAccession, Data_path):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{primaryAccession}-F1-model_v4.pdb"
    filename = f"AF-{primaryAccession}-F1-model_v4.pdb"
    filename = os.path.join(Data_path, filename)
    response = requests.get(url)
    if ("<Error><Code>NoSuchKey</Code>" in response.text):
       sys.exit("There does not exist a PDB file for that gene in AlphaFold")
    else:
       with open(filename, 'wb') as f: 
           f.write(response.content)
    return filename
