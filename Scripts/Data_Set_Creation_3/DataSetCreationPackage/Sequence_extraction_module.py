#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

from Bio import SeqIO
import numpy as np
import pandas as pd
import os

#####################################################################################################################################################################################################################

"""
Esta función recoge tanto el target obtenido para un gen de ensembl como el dataPath donde está almacenado el archivo. De esta manera, nos devolverá
el transcript que utilizan en la base de datos de CCDS que lo relacionará con el gen de interés. Esto es necesario dado que necesitamos ese id 
para buscar la secuencia de la proteína como la del ADN en las bases de datos de CCDS, dado que estas están en formato FASTA y su id corresponde
al CCDS_transcrip.
"""


def extract_CCDS(target, Data_Path):
    file = os.path.join(Data_Path, "CCDS2Sequence.current.txt")
    df   = pd.read_csv(file, sep = "\t", header = 0)
    ccds_transcript = df[df["nucleotide_ID"] == target]["#ccds"].values[0]
    return ccds_transcript

"""
Esta función recoge el nombre del archivo que contiene las secuencias, y el transcript a buscar. De esta manera, lo primero que hace es recoger
todos las secuencias exitentes en dicho archivo en una lista, y después analiza el nombre que se le ha dado a cada secuencia y lo compara con
el target introducido para obtener la secuencia de interés.
"""

def searching_sequence(query_sequence_file, target):
    query_sequence   = [i for i in SeqIO.parse(query_sequence_file, 'fasta')]
    sequence_interes = []
    for sequence in query_sequence:
        if (sequence.name.split("|")[0] == target):
            sequence_interes = sequence
            break
    return sequence_interes


"""
Esta función recoge las secuencias correspondientes a los nucleotidos y las proteínas gracias al módulo sequence_extr_mod, y las prepara para devolvernos un dataframe, 
relacionando los codones con el aminoácido correspondiente.
"""

def prepare_dataFrame(nucleotide_Sequence, protein_Sequence):
    nucleotide = [nucleotide_Sequence[i:i+3] for i in range(0, len(nucleotide_Sequence), 3)]
    protein    = [*protein_Sequence]
    while len(protein) != len(nucleotide):
        protein.append(np.nan)
    col = {"Nucleotide" : nucleotide, "Protein" : protein}
    df  = pd.DataFrame(data = col)
    return df
