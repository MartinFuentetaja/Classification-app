#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd
import numpy as np
import math
import subprocess, os
from Bio import AlignIO, SeqIO


#####################################################################################################################################################################################################################

"""
Esta funcíon recoge dos archivos de texto: in_file y out_file. in_file corresponde a un archivo tipo FASTA y es donde estarán almacenadas todas las secuencias de proteínas que 
queramos alinear. Por otro lado, out_file es el archivo donde se almacenarán las secuencias alineadas. De esta manera, utilizaremos el programa MUSCLE para construir el MSA.
Así mismo, esta función analiza el tamaño del archivo in_file, dado que si este excede los 40.0KB la manera de ejecutar el programa MUSCLE será diferente, dado que utilizaremos
"-super5" que en principio alinea con mayor rapidez un mayor número de secuencias. El output corresponde al out_file que contiene las secuencias alineadas.
"""

def obtain_MSA_MUSCLE(in_file, out_file):
    file_stats = os.stat(in_file)
    file_size  = file_stats.st_size * 0.008
    muscle_exe = "/home/einstein/martin/git_enviroment/database_toolkits/muscle5.1.linux"
    if (file_size < 40.0):
        muscle_result = subprocess.check_output([muscle_exe, "-align", in_file, "-output", out_file])
    else:
        muscle_result = subprocess.check_output([muscle_exe, "-super5", in_file, "-output", out_file])
    return out_file

"""
Esta función recoge un array y devuelve otro array el cual contiene la conservación del residuo por posición. Así pues, está función calcula la "Shanon_entropy", pero tiene en cuenta
las penalizaciones de los gaps ("-"). Para ello, en vez de calcular la probabilidad como en los otros casos (número de un aminoácido/número de secuencias), lo que se hace es 
dar al gap siempre la misma probabilidad (1 / número de secuencias) * el número de apariciones en dicha columna. 
"""

def Shanon_entropy(x):
    n = np.sum(x)
    probability = x[:-1] / n 
    prob_log    = np.ma.log(probability)
    probability = np.append(probability, np.array(1.0 * x[-1]/ n))
    prob_log    = np.append(prob_log, np.ma.log(np.array(1.0 / n)))
    entropy     = -np.dot(probability, prob_log)
    normalized_entropy = 1.0 - (entropy / math.log(n))
    return normalized_entropy

"""
Esta función recoge el aliniamiento del MSA y construye un DataFrame, donde cada fila corresponde a cada columna del aliniamiento, y hay 21 columnas, 20 corresponden a los aminácidos
y 1 al gap "-". De esta manera, realiza un histograma por cada columna, es decir, cuenta el número de apariciones de cada aminoácido y el de los gaps para poder hacer el cáclulo 
del "Shanon´s entropy" por cada columna del aliniamiento.
"""

def alnSiteCompositionDF(aln, characters="ACDEFGHIKLMNPQRSTVWY-"):
  alnRows = aln.get_alignment_length()
  compDict = {char:[0]*alnRows for char in characters}
  for record in aln:
    header = record.id
    seq = record.seq
    for aaPos in range(len(seq)):
      aa = seq[aaPos]
      if aa in characters:
        compDict[aa][aaPos] += 1    
  df = pd.DataFrame.from_dict(compDict)
  return df

"""
Esta función calcula la conservación del residuo por columnas del MSA, es decir, calcula cada conservación de cada columna y la añade a una lista. Una vez termine el análisis 
devuelve la lista con las conservaciones por cada posición. A pesar de esto, hay que tener en cuenta que la longitud de dicha lista puede no coincidir con la longitud de la secuencia
que estemos trabajando, como es el caso del KCNQ3, es decir, la secuencia es de 872 pero al construir el MSA esta cambia. Por lo tanto, hay que hacer ciertas modificaciones a esta.
Para ello utilizaremos la función adjusting_residue_conservation(). 
"""

def residue_conservation(in_file, out_file, returned_length):
    out_file = obtain_MSA_MUSCLE(in_file, out_file)
    align = AlignIO.read(out_file, "fasta")
    df    = alnSiteCompositionDF(align)
    #df.to_excel("b.xlsx", index = 0)
    residue_conserv = []
    rows, columns   = df.shape
    for i in range(rows):
        x = np.array(df.loc[i].tolist())
        residue_conserv.append(Shanon_entropy(x))
    residue_conserv = adjusting_residue_conservation(in_file, out_file, returned_length, residue_conserv)
    return residue_conserv

"""
Esta función recoge el in_file, out_file y returned_length. Esta última corresponde a la longitud de aminoácidos que tiene la secuencia que hemos obtenido de AlphaFold. De esta forma,
lo primero que encuentra es el id que tiene la misma longitud de aminoácidos que la de AlphaFold. Una vez lo tenemos, acudimos al out_file, leemos los aliniamientos, y extraemos la 
la secuencia correspondiente. Una vez logrado esto, encotramos las posiciones donde se están los gaps, y las extraemos de la lista del residue_conserv. De esta manera, la longitud
de la lista y de la secuencia serán la misma.
"""

def adjusting_residue_conservation(in_file, out_file, returned_length, residue_conserv):
    id = find_sequence_id(in_file, returned_length)
    align = AlignIO.read(out_file, "fasta")
    for al in align:
        if (al.id == id):
            seq = al.seq
            break
        else:
            continue
    seq = [*seq]
    ind = encontrar_posiciones(seq)
    residue_conserv = borrar_elementos(residue_conserv, ind)
    #os.remove(in_file)
    #os.remove(out_file)
    return residue_conserv

def find_sequence_id(in_file, wanted_length):
    align = SeqIO.parse(in_file, "fasta")
    for al in align:
        if (len(al.seq) == wanted_length):
            id = al.id
            break
        else:
            continue
    return id


def encontrar_posiciones(lista):
    return [i for i, e in enumerate(lista) if e == "-"]


def borrar_elementos(lista, posiciones):
    for indice in sorted(posiciones, reverse=True):
        del lista[indice]
        #lista[indice] = np.nan
    return lista
