#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import math

import DataSetCreationPackage.Amino_acid_calculation_module as amino_mod
import DataSetCreationPackage.clinvar_module                as cv_mod


#####################################################################################################################################################################################################################

"""
Esta función recoge tanto el dataFrame de AlphaFold como el de ClinVar, y devuelve una lista con el pLDDT. Lo que sabemos, es que el pLDDT está en la 
columna "b_factor" del DataFrame de AlphaFold. Por lo tanto, lo que hacemos es comparar la posición en la secuencia de la proteína que indica el cambio
de un amino ácido del DataFrame de ClinVar con el de AlphaFold. Si encontramos coincidencias entonces obtiene el pLDDT para dicha posición y lo almacena
en una lista. Sino, introduce np.nan, es decir, lo deja en blanco.
"""

def comparing_AlphaFold_Clinvar_pLDDT(df_AlphaFold, df_Clinvar):
    df_AlphaFold_unique = df_AlphaFold.drop_duplicates(subset = ['residue_number'], keep ='first')
    merged_df = pd.merge(df_Clinvar, df_AlphaFold_unique, left_on='position', right_on='residue_number', how='left')
    pLDDT_list = merged_df['b_factor'].tolist()
    pLDDT_list = [np.nan if pd.isnull(plddt) else plddt for plddt in pLDDT_list]
    return pLDDT_list


"""
Esta función recoge la lista dvuelta por la función "getting_Secondary_Structure()" que contiene la estrcutura secundaria de la proteína, y el DataFrame
de ClinVar. Lo que sabemos es que cada índice de la lista de "secondary_structure_list" corresponde a una posición en la secuencia de la proteína:
índcie      Position
  0             1
  1             2
  2             3

De esta manera, lo unico que devemos hacer es leer la columna "position" del DataFrame de ClinVar, y buscar el indice correspondiente en la lista de
la estructura secundaria. Una vez hemos hecho esto, el término correspondiente a dicho índice se almacena en una lista. A pesar de esto, se hace un 
control de longitudes por si se diera el caso que la posición buscada en el ClinVar este fuera del rango de la lista. De ser así, se añadirá un 
vacío (np.nan). Tambíen se utiliza para residue_conserv
"""

def comparing_Clinvar_obtain_list(compared_list, df_Clinvar):
    df_positions = df_Clinvar["position"].astype('int')
    list_mod = df_positions.apply(lambda df: compared_list[df-1] if df <= len(compared_list) else np.nan).tolist()
    return list_mod


"""
Esta función recoge el DataFrame del Clinvar y el de amino ácido (que es fijo). De esta manera, utilizando el módulo de "Amino_acid_calculation_module"
obtiene los descriptores correspondientes y los almacena en listas para añadirlos posteriormente al dataFrame.
"""

def adding_amino_acid_features(df, df_amino):
    rows, columns = df.shape
    charge_list, polarity_list, aroma_list, d_size_list, d_vol_list, d_pol_e_list, d_ip_e_list, d_hf_e_list, d_msa_list = [], [], [], [], [], [], [], [], []
    for i in range(rows):
        charge_change, polarity_change, aroma_change, d_size, d_vol, d_pol_e, d_ip_e, d_hf_e, d_msa = amino_mod.amino_acid_comparision(df.loc[i, "initial amino acid"], df.loc[i, "final amino acid"], df_amino)
        charge_list.append(charge_change)
        polarity_list.append(polarity_change)
        aroma_list.append(aroma_change)
        d_size_list.append(d_size)
        d_vol_list.append(d_vol)
        d_pol_e_list.append(d_pol_e)
        d_ip_e_list.append(d_ip_e)
        d_hf_e_list.append(d_hf_e)
        d_msa_list.append(d_msa)
    return charge_list, polarity_list, aroma_list, d_size_list, d_vol_list, d_pol_e_list, d_ip_e_list, d_hf_e_list, d_msa_list

"""
Esta función recoge el DataFrame con el fin de obtener sus dimensiones, y devuelve las tres listas vacías. Estas se añadirán posteriormente al excel.
Esto está de esta forma dado que son descriptores que Alvaro los define a mano. 
"""
def adding_Structural_features(df):
    rows, columns = df.shape
    domain_list   = [np.nan]*rows
    function_list = [np.nan]*rows
    region_list   = [np.nan]*rows
    return domain_list, function_list, region_list

"""
Esta función recoge tanto el DataFrame de ClinVar como el de Lovd3. Una vez leídos, comparará la información relevante
a la clasificación clínica por si hay discrepancias. De esta manera, escogerá los elementos de la columna Protein y DNA change 
que sean iguales, y dividirá el dataFrame principal en dos: coincidente y no. Por otro lado, si ha encontrado coincidencias analizará las posibles
diferencias que pueda haber a la hora de haber clasificado la condición clínica. Por lo tanto, se da más relevancia a la 
clasificación donde no se muestre la palabra "likely". En caso de estar en las dos no se realizará ninguna modificación. Finalmente, llama a la función
Lovd3_not_coincidences para añadir la info no coincidente a la coincidente.
"""

def compare_clinvar_lovd3(df_clinvar, df_lovd3, length):
    df_lovd3["DNA change (cDNA)"]        = df_lovd3["DNA change (cDNA)"].str.strip()
    df_lovd3["Protein"]                  = df_lovd3["Protein"].str.strip()
    df_clinvar["DNA change"]             = df_clinvar["DNA change"].str.strip()
    df_clinvar["Protein change (compl)"] = df_clinvar["Protein change (compl)"].str.strip()

    df = pd.merge(df_clinvar, df_lovd3, left_on = ["Protein change (compl)", "DNA change"], right_on = ["Protein", "DNA change (cDNA)"], how = "outer")
    
    df_coincidences     = df[~df["Name"].isnull().values].reset_index(drop = True)
    df_not_coincidences = df[df["Name"].isnull().values & ~df["Clinical classification"].isnull().values].reset_index(drop = True)
    df_coincidences["Clinical classification"]  = df_coincidences["Clinical classification"].fillna("PASS")
    
    condition = (~df_coincidences["Clinical classification"].str.contains("PASS") &
    (
        (
            (df_coincidences["Clinical significance (Last reviewed)"].str.lower().str.contains("likely")) &
            (~df_coincidences["Clinical classification"].str.lower().str.contains("likely"))
            ) | 
        (df_coincidences["Clinical significance (Last reviewed)"].str.lower().str.contains("uncertain|not provided|no interpretations", case=False))
            )
                 )
    df_coincidences.loc[condition, "Clinical significance (Last reviewed)"] = df_coincidences.loc[condition, "Clinical classification"]
    df_coincidences = Lovd3_not_coincidences(df_coincidences, df_not_coincidences, length)
    return df_coincidences

"""
Esta función es complementaria de la anterior. De esta forma, recoge los dataFrames coincidentes y no, y une la info de la no coincidente a la coincidente.
"""

def Lovd3_not_coincidences(df_coincidences, df_not_coincidences, length):
    rows, columns                        = df_not_coincidences.shape
    variation_name_ref                   = df_coincidences.loc[0, "Name"]
    DNA_change_ref                       = df_coincidences.loc[0, "DNA change"]
    Protein_change_ref                   = df_coincidences.loc[0, "Protein change (compl)"]
    for i in range(rows):
        DNA_change     = df_not_coincidences.loc[i, "DNA change (cDNA)"]
        Protein_change = df_not_coincidences.loc[i, "Protein"]
        variation_name = variation_name_ref.replace(DNA_change_ref, DNA_change).replace(Protein_change_ref, Protein_change)
        position, init_amino_mod, fin_amino_mod, change, compl_change = cv_mod.extracting_protein_change(Protein_change)
        nucl_change, nucl_pos, initial_nucl, final_nucl               = cv_mod.obtain_nucleotido_change(DNA_change)
        if ((init_amino_mod != "*" and fin_amino_mod != "*") and int(position) <= int(length)):
            variation_complete = {"Name" : variation_name, "Gene(s)" : df_coincidences.loc[0, "Gene(s)"], "DNA change" : DNA_change, "DNA position" : nucl_pos, "initial nucleotide" : initial_nucl,
                                    "final nucleotide" : final_nucl, "Protein change (compl)" : Protein_change, "Protein change" : change, "position" : position, "initial amino acid" : init_amino_mod, 
                                    "final amino acid" : fin_amino_mod, "Condition(s)" : np.nan, "Clinical significance (Last reviewed)" : df_not_coincidences.loc[i, "Clinical classification"], "Review status" : np.nan}
            df_coincidences.loc[len(df_coincidences)] = variation_complete
    df_coincidences = df_coincidences.drop(columns = ["DNA change (cDNA)", "Protein", "Clinical classification"])
    df_coincidences = df_coincidences.sort_values(by = ["position"], ascending = True, na_position = "last").reset_index(drop=True)
    return df_coincidences


"""
Esta función recoge el dataFrame resultante de GenomAD_and_ClinVar_comparision_1, y realizará comparaciones entre las columnas clínicas de la existente con la procedente de GenomAD. 
Después de realizar comparaciones se quedará con el mejor caso.
"""

def GenomAD_and_ClinVar_coincidences_comparision(df_coincidences):
    df_coincidences["ClinVar Clinical Significance"] = df_coincidences["ClinVar Clinical Significance"].fillna("PASS")
    conditions = (~df_coincidences["ClinVar Clinical Significance"].str.contains("PASS") & df_coincidences["ClinVar Clinical Significance"].str.lower().str.contains("uncertain|not provided"))
    df_coincidences.loc[conditions, "ClinVar Clinical Significance"] = "Benign"
    conditions_general = (
        ~df_coincidences["ClinVar Clinical Significance"].str.contains("PASS") & 
            (
                ((df_coincidences["Clinical significance (Last reviewed)"].str.lower().str.contains("uncertain|not provided|no interpretation")) & (~df_coincidences["ClinVar Clinical Significance"].str.lower().str.contains("uncertain|conflicting|not provided"))) |
                ((df_coincidences["Clinical significance (Last reviewed)"].str.lower().str.contains("likely")) & (~df_coincidences["ClinVar Clinical Significance"].str.lower().str.contains("likely|uncertain|conflicting|not provided"))) |
                (~(df_coincidences["Clinical significance (Last reviewed)"].str.lower().str.contains("pathogenic")) & (df_coincidences["ClinVar Clinical Significance"].str.lower().str.contains("benign"))) |
                (~(df_coincidences["Clinical significance (Last reviewed)"].str.lower().str.contains("benign")) & (df_coincidences["ClinVar Clinical Significance"].str.lower().str.contains("pathogenic")))
                )
        )
    df_coincidences.loc[conditions_general, "Clinical significance (Last reviewed)"] = df_coincidences.loc[conditions_general, "ClinVar Clinical Significance"]
    df_coincidences = df_coincidences.drop(columns = ["ClinVar Clinical Significance", "Protein Consequence", "Transcript Consequence"]).reset_index(drop = True)
    return df_coincidences
    

"""
Esta función recoge el dataFrame resultante de la combinación del existente con el obtenido de GenomAD. Extraerá dos dataFrames, uno correspondiente a las coincidencias obtenidas 
y el otro a las no coincidentes. Esto es necesario para utilizarlo en la función "GenomAD_and_ClinVar_comparision_1". 
"""

def extact_non_null_elements(df):
    df_coincidences     = df[~df['Name'].isnull().values].reset_index(drop=True)
    df_not_coincidences = df[df['Name'].isnull().values].reset_index(drop=True)
    df_not_coincidences.loc[(df_not_coincidences["ClinVar Clinical Significance"].isnull() | df_not_coincidences["ClinVar Clinical Significance"].str.lower().str.contains("uncertain|not provided")), "ClinVar Clinical Significance"] = "Benign"
    return df_coincidences, df_not_coincidences

"""
Esta función recogerá el dataFrame resultante de la unión de ClinVar y Lovd3, y lo unirá con el de GenomAD. De esta manera, dividirá los casos coincidente y no,
y utilizará una función externa para hacer una comparación en la clasificación clínica para los casos coincidentes. Para los no coincidentes, les asignará a todos
como benignos y los áñadirá al coincidente.
"""

def GenomAD_and_ClinVar_comparision(df_final, df_GenomAD, length):
    df_GenomAD["Transcript Consequence"] = df_GenomAD["Transcript Consequence"].str.strip()
    df_GenomAD["Protein Consequence"]    = df_GenomAD["Protein Consequence"].str.strip()
    df = pd.merge(df_final, df_GenomAD, left_on = ["Protein change (compl)", "DNA change"], right_on = ["Protein Consequence", "Transcript Consequence"], how = "outer")
    df_coincidences, df_not_coincidences = extact_non_null_elements(df)
    df_coincidences = GenomAD_and_ClinVar_coincidences_comparision(df_coincidences)
    rows, columns                        = df_not_coincidences.shape
    variation_name_ref                   = df_coincidences.loc[0, "Name"]
    DNA_change_ref                       = df_coincidences.loc[0, "DNA change"]
    Protein_change_ref                   = df_coincidences.loc[0, "Protein change (compl)"]
    for i in range(rows):
        DNA_change     = df_not_coincidences.loc[i, "Transcript Consequence"]
        Protein_change = df_not_coincidences.loc[i, "Protein Consequence"]
        variation_name = variation_name_ref.replace(DNA_change_ref, DNA_change).replace(Protein_change_ref, Protein_change)
        position, init_amino_mod, fin_amino_mod, change, compl_change = cv_mod.extracting_protein_change(Protein_change)
        nucl_change, nucl_pos, initial_nucl, final_nucl = cv_mod.obtain_nucleotido_change(DNA_change)
        if ((init_amino_mod != "*" and fin_amino_mod != "*") and int(position) <= int(length)):
            variation_complete = {"Name" : variation_name, "Gene(s)" : df_coincidences.loc[0, "Gene(s)"], "DNA change" : DNA_change, "DNA position" : nucl_pos, "initial nucleotide" : initial_nucl,
                                    "final nucleotide" : final_nucl, "Protein change (compl)" : Protein_change, "Protein change" : change, "position" : position, "initial amino acid" : init_amino_mod, 
                                    "final amino acid" : fin_amino_mod, "Condition(s)" : np.nan, "Clinical significance (Last reviewed)" : df_not_coincidences.loc[i, "ClinVar Clinical Significance"], "Review status" : np.nan}
            df_coincidences.loc[len(df_coincidences)] = variation_complete
    df_coincidences = df_coincidences.sort_values(by = ["position"], ascending = True, na_position = "last").reset_index(drop=True)
    return df_coincidences
