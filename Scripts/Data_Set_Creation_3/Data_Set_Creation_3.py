#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd
import numpy as np
import math
import os, shutil, sys
import mdtraj as md
import requests
import time
from Bio import AlignIO, SeqIO
from biopandas.pdb import PandasPdb

import DataSetCreationPackage.clinvar_module                as clv_mod
import DataSetCreationPackage.LOVD3_file_module             as lovd3_file_mod
import DataSetCreationPackage.GenomAD_file_module           as GenomAD_file_mod
import DataSetCreationPackage.Amino_acid_calculation_module as amino_mod
import DataSetCreationPackage.Calculations_df_frames_mod    as calc_df_mod
import DataSetCreationPackage.One_hot_module                as one_hot_mod
import DataSetCreationPackage.residue_conserv_module        as res_conserv_mod
import DataSetCreationPackage.MSA_HHBLITS 		            as msa_hhblits
import DataSetCreationPackage.File_obtention_module_2_0     as file_obtention_mod
import DataSetCreationPackage.clinvar_database_module       as ClinVar_data_mod
import DataSetCreationPackage.WebPackage.web_Ensembl_module as web_ensembl_mod
import DataSetCreationPackage.Sequence_extraction_module    as sequence_extr_mod
import DataSetCreationPackage.MTR_calc_module               as MTR_calc_mod
import DataSetCreationPackage.argparse_decision             as argparse_decision


#####################################################################################################################################################################################################################

"""
Esta función recoge el nombre de un archivo y devuelve el path donde se encuentra.
"""

def controlling_if_file_exists(results_file):
    pwd  = os.path.dirname(__file__)
    path = os.path.join(pwd, results_file)
    return os.path.exists(path)

"""
Esta función recoge el nombre de un directorio que queremos crear. La función evalua si el directorio existe o no, y de no existir, lo crea.
"""

def directory_creation(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

"""
Esta función recoge el archivo de AlphaFold (.pdb) y el status (True/False) que indica si queremos una descripción con mayor detalle. De esta manera,
lee el archivo de AlphaFold, y directamente extrae la estructura secundaria. El resultado nos lo devuelve en un array ([[]]).

If Status == False

The DSSP assignment codes are:

-‘H’ : Alpha helix

-‘B’ : Residue in isolated beta-bridge

-‘E’ : Extended strand, participates in beta ladder

-‘G’ : 3-helix (3/10 helix)

-‘I’ : 5 helix (pi helix)

-‘T’ : hydrogen bonded turn

-‘S’ : bend

-‘ ‘ : Loops and irregular elements

If status == True

The simplified DSSP codes are:

-‘H’ : Helix. Either of the ‘H’, ‘G’, or ‘I’ codes.

-‘E’ : Strand. Either of the ‘E’, or ‘B’ codes.

-‘C’ : Coil. Either of the ‘T’, ‘S’ or ‘ ‘ codes.

"""

def getting_Secondary_Structure(AlphaFold_file, status):
    pdb = md.load_frame(AlphaFold_file, 0)
    secondary_struct = md.compute_dssp(pdb, simplified=status)
    return secondary_struct

"""
Esta función calcula la conservación del residuo por columnas del MSA, es decir, calcula cada conservación de cada columna y la añade a una lista. Una vez termine el análisis
devuelve la lista con las conservaciones por cada posición. A pesar de esto, hay que tener en cuenta que la longitud de dicha lista puede no coincidir con la longitud de la secuencia
que estemos trabajando, como es el caso del KCNQ3, es decir, la secuencia es de 872 pero al construir el MSA esta cambia. Por lo tanto, hay que hacer ciertas modificaciones a esta.
Para ello utilizaremos la función adjusting_residue_conservation().
"""

def residue_conservation(in_file, out_file, cpu_number, Save_path, UniRefPath, BfdPath):
    if len(out_file.split(".")) < 2:
       out_file = msa_hhblits.obtain_MSA_HHBLITS(in_file, out_file, cpu_number, Save_path, UniRefPath, BfdPath)
    align = AlignIO.read(out_file, "fasta")
    df    = msa_hhblits.alnSiteCompositionDF(align)
    #df.to_excel("b.xlsx", index = 0)
    residue_conserv = []
    rows, columns   = df.shape
    for i in range(rows):
        x = np.array(df.loc[i].tolist())
        residue_conserv.append(msa_hhblits.Shanon_entropy(x))
    return residue_conserv

"""
Esta función recoge la base de datos construída y realiza una codificación de tipo one_hot, es decir para ciertos valores se les asigna un 1 mientras que para otros un 0. Las funciones utilizadas se encuentran en el módulo one_hot_mod.
"""

def one_hot_file(filename):
    data                     = pd.read_excel(filename, header = 0)
    amino_encoded            = one_hot_mod.one_hot_amino(data)
    second_structure_encoded = one_hot_mod.one_hot_second_structure(amino_encoded)
    amino_changes_encoded    = one_hot_mod.one_hot_general(second_structure_encoded, ["Charge_change", "Polarity_change", "Aromaticity_change"])
    my_label_encoded         = one_hot_mod.my_label_modification(amino_changes_encoded)
    df_final                 = one_hot_mod.data_Frame_modification(my_label_encoded)
    filename_hot             = os.path.basename(filename).split(".")[0] + "_encoded.xlsx"
    dirname                  = os.path.dirname(filename)
    df_final.to_excel(os.path.join(dirname, filename_hot), index = False)
    return filename_hot

"""
Esta función recoge el nombre del gen y el PATH donde se almacena toda la Data. De esta manera el primer paso es obtener el target que relaciona
el nombre del gen con el id proporcionado por ensembl. Una vez obtenido esto se obtiene el transcript que relaciona el id de ensembl con el de CCDS.
Finalmente, sabiendo cual es el id de CCDS se extraen las secuencias correspondientes a ese gen, y finalmente devuelve un dataframe que relaciona
los codones con sus correspondientes amino ácidos.
"""

def codon_protein_df_initializing(gene_name, lenght,  Data_path):
    target                 = web_ensembl_mod.extract_id_ensembl(gene_name, lenght)
    ccds_transcript        = sequence_extr_mod.extract_CCDS(target, Data_path)
    nucleotide_file        = os.path.join(Data_path, "CCDS_nucleotide.current.fna")
    protein_file           = os.path.join(Data_path, "CCDS_protein.current.faa")
    nucleotide_sequence    = str(sequence_extr_mod.searching_sequence(nucleotide_file, ccds_transcript).seq)
    protein_sequence       = str(sequence_extr_mod.searching_sequence(protein_file, ccds_transcript).seq)
    df_pro_nuc             = sequence_extr_mod.prepare_dataFrame(nucleotide_sequence, protein_sequence)
    df_pro_nuc["position"] = df_pro_nuc.index + 1
    return df_pro_nuc


"""
Esta función recoge dos dataFrames: df_CodonProtein y df_GenomAD, donde el primero contiene la relación Codón Proteína, y el segundo es la unión
entre la info obtenida de GenomAD_V2 y GenomAD_V3_non_V2. Además, recoge también el PATH donde se almacena el archivo CodigoGenetico para poder
extraer el número tanto de sinónimos y MissNotpos por cada Codón y Amino ácido. De esta manera, lo primero que se realiza es una unión entre 
df_CodonProtein y df_GenomAD por posición, para relacionar cada codón con posibles variaciones, y sus Allele Numbers y Count correspondientes.
Una vez realizado esto, se hace una segunda unión entre la base de datos del código Genético y el dataFrame resultante de la unión anterior. 
De esta forma, relacionamos los valores de Missense y Synonimous con sus correspondientes codones. Una vez recogida toda la info somos capaces
de extraer el MTI, pero primero se calcula la frecuencia, después el MTR, y finalmente el MTI.
"""

def MTI_calculation(df_CodonProtein, df_GenomAD, Data_path):
    length           = df_CodonProtein.shape[0]
    GeneticCode_file = os.path.join(Data_path, "CodigoGenetico.xlsx")
    df_merge         = pd.merge(df_CodonProtein, df_GenomAD, on = ["position"], how = "left")
    df_merge         = Genetic_code_df(GeneticCode_file, df_merge)
    MTR              = MTR_calc_mod.Frequency_calculation(df_merge, length)
    MTI              = MTR_calc_mod.MTI_ponderado(MTR, 3)
    return MTI

"""
Esta función recoge el dataFrame creado, y compara las columnas de "Codigo" del CodigoGentico.xlsx con"Codons" del dataFrame para encontrar coincidencias, 
y de esta manera añadir los sinónimos correspondientes.
"""

def Genetic_code_df(GeneticFile, df_base):
    df_genetic = pd.read_excel(GeneticFile, header = 0)
    df = pd.merge(df_base, df_genetic[["Codigo", "LetraAA", "MissensePos", "SynonymousPos", "MissNoStopPos"]], left_on = "Nucleotide", right_on = "Codigo", how = "left")
    df = df.drop(columns = ["Codigo", "LetraAA"]).reset_index(drop=True)
    return df


"""
Esta función recoge los ficheros obtenidos de GenomAD, e inicializa cada DataFrame. Una vez realizado esto se hace una unión por posición, proteína
y DNA, para relacionar todos los casos, y poder saber el valor correspondiente a cada posición de Alle Count y NUmber por cada fichero.
"""

def Initializing_GenomAD(GenomAD_V2_file, GenomAD_V3_file):
    df_V2 = GenomAD_file_mod.GenomAD_modification(GenomAD_V2_file)
    df_V2 = GenomAD_file_mod.GenomAD_second_modification(df_V2)
    df_V3 = GenomAD_file_mod.GenomAD_modification(GenomAD_V3_file)
    df_V3 = GenomAD_file_mod.GenomAD_second_modification(df_V3)
    df_GenomAD = pd.merge(df_V2[["position", "Protein Consequence", "Transcript Consequence", "Allele Number", "Allele Count"]], df_V3[["position", "Protein Consequence", "Transcript Consequence", "Allele Number", "Allele Count"]], on = ["position","Protein Consequence", "Transcript Consequence"], how = "outer", sort=True)
    return df_GenomAD

def Initializing_GenomAD_2(GenomAD_V2_control_file, GenomAD_V2_non_neuro_file):
    df_V2_control   = GenomAD_file_mod.GenomAD_modification(GenomAD_V2_control_file)
    df_V2_control   = GenomAD_file_mod.GenomAD_second_modification(df_V2_control)
    df_V2_non_neuro = GenomAD_file_mod.GenomAD_modification(GenomAD_V2_non_neuro_file)
    df_V2_non_neuro = GenomAD_file_mod.GenomAD_second_modification(df_V2_non_neuro)
    df_GenomAD = pd.merge(df_V2_control[["Protein Consequence", "Transcript Consequence", "ClinVar Clinical Significance"]], df_V2_non_neuro[["Protein Consequence", "Transcript Consequence", "ClinVar Clinical Significance"]], on = ["Protein Consequence", "Transcript Consequence", "ClinVar Clinical Significance"], how = "outer", sort=True)
    return df_GenomAD

"""
Esta función recoge el nombre del arcvio config.txt que contiene el PATH donde están almacenadas todas las bses de datos. De esta manera,
el usuario escribirá en dicho archivo donde contiene el fichero DATA, y el programa lo leerá. Además, tiene dos controles por si el PATH escrito
no existe, o por si no hay nada escrito.
"""

def data_path(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        try:
            for line in lines:
                if (line.strip().split("=")[0] == ">DATA PATH"):
                    path = line.strip().split("=")[1]
                    if (os.path.exists(path)):
                        DataPath = path
                    else:
                        sys.exit("The path provided for Data in config.txt does not exist.")
                elif (line.strip().split("=")[0] == ">UNIREF PATH"):
                    path = line.strip().split("=")[1]
                    if (os.path.exists(os.path.dirname(path))):
                        UniRefPath = path
                    else:
                        sys.exit("The path provided for UniRef in config.txt does not exist.")
                elif (line.strip().split("=")[0] == ">BFD PATH"):
                    path = line.strip().split("=")[1]
                    if (os.path.exists(os.path.dirname(path))):
                        BfdPath = path
                    else:
                        sys.exit("The path provided for Bfd in config.txt does not exist.")
            return DataPath, UniRefPath, BfdPath
        except UnboundLocalError:
            sys.exit("The congif.txt file is badly configurated. Data path must be write It as: >DATA PATH=PATH... & >UNIREF PATH=... & >BFD PATH=...")



#Introducing inputs

args = argparse_decision.parser_argument()

#Inicializando los ficheros
Data_path, UniRefPath, BfdPath = data_path(args.config_path)
Current_path                   = os.path.dirname(Data_path)

#Creando el fichero de los resultados.

if args.filename:
    ClinVar_file              = args.filename[1]
    AlphaFold_file            = args.filename[0]
    in_file                   = args.filename[2]
    Lovd3_file                = args.filename[3]
    GenomAD_V2_file           = args.filename[4]
    GenomAD_V2_control_file   = args.filename[5]
    GenomAD_V2_non_neuro_file = args.filename[6]
    GenomAD_V3_file           = args.filename[7]
    gene_name                 = in_file.split(".")[0]
else:
    gene_name   = args.gen_name
    length      = args.length

    #Downloading files:
    #ClinVar_data_mod.date_control()
    if (args.web_decision == "Internet" or args.web_decision == "Int" or args.web_decision == "I"):
       AlphaFold_file, returned_length, in_file, ClinVar_file, lovd3_file, GenomAD_V2_file, GenomAD_V2_control_file, GenomAD_V2_non_neuro_file, GenomAD_V3_file = file_obtention_mod.obtaining_files_web(gene_name, length, Data_path)
    else:
       AlphaFold_file, returned_length, in_file, ClinVar_file, lovd3_file, GenomAD_V2_file, GenomAD_V2_control_file, GenomAD_V2_non_neuro_file, GenomAD_V3_file = file_obtention_mod.obtaining_files_web_DataBase(gene_name, length, Data_path)



#ClinVar modifications:
df       = clv_mod.NM_modification_clinvar(ClinVar_file)
df_mod   = clv_mod.NM_second_modification_clinvar(df)
df_final = clv_mod.Protein_modification_clinvar_2(df_mod)

#Lovd3 modifications:
df_lovd3 = lovd3_file_mod.modify_Lovd3(lovd3_file)
df_lovd3 = lovd3_file_mod.second_modification_Lovd3(df_lovd3)

#GenomAD control and non_neuro initializing:
df_GenomAD_control_non_neu = Initializing_GenomAD_2(GenomAD_V2_control_file, GenomAD_V2_non_neuro_file)

#Combining Lovd3 and Clinvar dataframes:
df_final = calc_df_mod.compare_clinvar_lovd3(df_final, df_lovd3, length)

#Combining with GenomAD
df_final = calc_df_mod.GenomAD_and_ClinVar_comparision(df_final, df_GenomAD_control_non_neu, length)

rows, columns = df_final.shape

#Una vez realizadas las modificaciones creamos el dataFrame para AlphaFold
pdb_df        = PandasPdb().read_pdb(AlphaFold_file)
atom_df       = pdb_df.df["ATOM"]

#Ahora podemos obtener el pLDDT y el secondary_Structure
pLDDT_list    = calc_df_mod.comparing_AlphaFold_Clinvar_pLDDT(atom_df, df_final)
if (args.status == "Extend"):
    secondary_structure_list = getting_Secondary_Structure(AlphaFold_file, False)[0]
else:
    secondary_structure_list = getting_Secondary_Structure(AlphaFold_file, True)[0]
try: 
    Second_Structure_list    = calc_df_mod.comparing_Clinvar_obtain_list(secondary_structure_list, df_final)
except IndexError:
    print("")
    sys.exit("There has been an indexing problem. It might be related with the shape of the file. Please check If they chosen ones are the correct ones.")


#Calculando el residue conservation
out_file 	    = gene_name
directory_MSA       = directory_creation(os.path.join(Current_path, "msas"))
directory_MSA_gen   = directory_creation(os.path.join(directory_MSA, f"{gene_name}"))
if not os.path.exists(os.path.join(directory_MSA_gen, f"{gene_name}.aligned.gap.filtered.fas")):
   residue_conserv  = residue_conservation(in_file, out_file, args.cpu_number, directory_MSA_gen, UniRefPath, BfdPath)
else:
   out_file         = os.path.join(directory_MSA_gen, f"{gene_name}.aligned.gap.filtered.fas")
   residue_conserv  = residue_conservation(in_file, out_file, args.cpu_number, directory_MSA_gen, UniRefPath, BfdPath)
residue_conserv_mod = calc_df_mod.comparing_Clinvar_obtain_list(residue_conserv, df_final)

#Inicializando el data Frame para los aminoAcidos
amino_acid_file = os.path.join(Data_path, "Amino_Acid_Features.xlsx")
df_amino        = pd.read_excel(amino_acid_file, header = 0, index_col = 0)

"""
#Inicializando el dataFrame para las relaciones entre nucleotidos y proteinas
df_pro_nuc          = codon_protein_df_initializing(gene_name, Data_path)
codon_list          = calc_df_mod.comparing_Clinvar_obtain_list(df_pro_nuc["Nucleotide"].tolist(), df_final)
df_final.insert(loc = df_final.columns.get_loc("Allele Number_y") + 1, column = "Codons", value = codon_list)

#Añadiendo los sinónimos, MissensePos y MissNotStop
GeneticCode_file = os.path.join(Data_path, "CodigoGenetico.xlsx")
df_final         = Genetic_code_df(GeneticCode_file, df_final)
"""

#Obteniendo varios descriptores para los aminoAcidos
charge_list, polarity_list, aroma_list, d_size_list, d_vol_list, d_pol_e_list, d_ip_e_list, d_hf_e_list, d_msa_list = calc_df_mod.adding_amino_acid_features(df_final, df_amino)

#Structural Features
domain_list, function_list, region_list = calc_df_mod.adding_Structural_features(df_final)

#MTI calculation:
df_codon                = codon_protein_df_initializing(gene_name, length, Data_path)
codon_list              = calc_df_mod.comparing_Clinvar_obtain_list(df_codon["Nucleotide"].tolist(), df_final)
df_GenomAD_V2_V3_non_V2 = Initializing_GenomAD(GenomAD_V2_file, GenomAD_V3_file)
MTI                     = MTI_calculation(df_codon, df_GenomAD_V2_V3_non_V2, Data_path)
MTI_mod                 = calc_df_mod.comparing_Clinvar_obtain_list(list(MTI), df_final)

#Una vez realizado todos los cálculos añadimos los descriptores a la base de datos:

features = [pLDDT_list, Second_Structure_list, charge_list, polarity_list, aroma_list, d_size_list, d_vol_list, d_pol_e_list, d_ip_e_list,
            d_hf_e_list, d_msa_list, residue_conserv_mod, domain_list, function_list, region_list, MTI_mod, codon_list]

features_names = ["pLDDT", "Second_Structure", "Charge_change", "Polarity_change", "Aromaticity_change", "d_size", "d_vol", "d_pol_e",
                  "d_ip_e", "d_hf_e", "d_msa", "residue_conserv", "Domain", "Function", "Region", "MTI_31", "Codons"]

final_amino_loc = df_final.columns.get_loc("final amino acid")
for i in range(len(features)):
    df_final.insert(loc = final_amino_loc + 1 + i , column = features_names[i], value = features[i])

print("")
print("SUCCESS!")
print("")


df_final.to_excel(os.path.join(Current_path,f"{args.outfile}.xlsx"), index = False)

directory     = directory_creation(os.path.join(Current_path, "RESULTS"))
directory_gen = directory_creation(os.path.join(directory, f"{gene_name}"))
path          = os.path.join(Current_path, f"{args.outfile}.xlsx")
path_2        = os.path.join(directory_gen, f"{args.outfile}.xlsx")

bool  = True
while (bool == True):
    if not os.path.exists(path_2):
        shutil.move(path, path_2)
        bool = False
    else:
        os.remove(path_2)
        #path_2  = os.path.join(directory_gen, f"{args.outfile}_{str(count)}.xlsx")

#One-Hot:

One_hot_data   = one_hot_file(path_2)


        
print("")
print(f"Your database and encoded one have been created at {directory_gen}")
print("")

os.remove(AlphaFold_file)
os.remove(ClinVar_file)
os.remove(in_file)
os.remove(lovd3_file)
os.remove(GenomAD_V2_file)
os.remove(GenomAD_V2_control_file)
os.remove(GenomAD_V2_non_neuro_file)
os.remove(GenomAD_V3_file)
