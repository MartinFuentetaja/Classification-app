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



def obtain_MSA_HHBLITS(in_file, out_file, cpu_number, Save_path, UniRefPath, BfdPath):
    hhblits_exe         = "hhblits"
    hhfilter_exe        = "hhfilter"
    reformat            = "reformat.pl"
    out_path_a3m        = os.path.join(Save_path, f"{out_file}.a3m")
    out_path_hhr        = os.path.join(Save_path, f"{out_file}.hhr")
    out_filt_path       = os.path.join(Save_path, f"{out_file}.fil.a3m")
    out_align_filt_path = os.path.join(Save_path, f"{out_file}.fil.fas")
    
    subprocess.run([hhblits_exe, "-i", in_file, "-cpu", cpu_number, "-o", out_path_hhr, "-oa3m", out_path_a3m, "-n", "3", "-e", "0.001", "-maxseq", "1000000", "-realign_max", "100000", "-maxfilt", "100000", "-min_prefilter_hits", "1000", "-d", BfdPath, "-d", UniRefPath])
    subprocess.run([hhfilter_exe, "-i", out_path_a3m, "-o", out_filt_path])
    subprocess.run([reformat, "-r", out_filt_path, out_align_filt_path])
    filtered_file = filter_sequences_by_gap_percentage(out_align_filt_path, os.path.join(Save_path,f"{out_file}.aligned.gap.filtered.fas"), gap_threshold = 50)
    return filtered_file

#hhblits -i ./data/KCNQ3.fas -cpu 4 -oa3m /tmp/tmporap5gyi/output.a3m -o /dev/null -n 3 -e 0.001 -maxseq 1000000 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -d /bigdisk/scicomp/AlphaFold/DATA/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt -d /bigdisk/scicomp/AlphaFold/DATA/uniclust30/uniclust30_2018_08/uniclust30_2018_08


def filter_sequences_by_gap_percentage(msa_file, filtered_file = "aligned.gap.filtered.fas", gap_threshold = 50):
    sequences = []
    with open(msa_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = str(record.seq)
            gap_percentage = (sequence.count("-") / len(sequence)) * 100
            if gap_percentage <= gap_threshold:
                sequences.append(record)

    with open(filtered_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    return os.path.abspath(filtered_file)


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


