#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd

import DataSetCreationPackage.clinvar_module as clv_mod
#####################################################################################################################################################################################################################

"""
Esta función recogerá el nombre del archivo para el cual se creará un dataFrame. Para ello solo se seleccionarán las columnas de interes, y acto seguido se realiza una criva donde solo se escogen los casos donde la columna "VEP Annotation" es igual a "missense_variant".
"""

def GenomAD_modification(filename):
    df = pd.read_csv(filename, sep = ",",header = 0, usecols = ["Transcript", "Protein Consequence", "Transcript Consequence", "VEP Annotation", "ClinVar Clinical Significance", "Allele Count", "Allele Number"])
    df = df[df["VEP Annotation"] == "missense_variant"].reset_index(drop=True)
    return df

"""
Una vez obtenido el dataFrame esta función añadirá la columna position con la finalidad de ordenar porposición el dataFrame. 
"""

def GenomAD_second_modification(df):
    rows, columns = df.shape
    position_list = []
    for i in range(rows):
        position, init_amino_mod, fin_amino_mod, change, compl_change = clv_mod.extracting_protein_change(df.loc[i, "Protein Consequence"])
        if (fin_amino_mod == "*" or init_amino_mod == "*"):
            df = df.drop(i)
        else:
            position_list.append(position)
    Protein_change_loc = df.columns.get_loc("Protein Consequence")
    df.insert(loc = Protein_change_loc + 1, column = "position", value = position_list)
    df = df.sort_values(by = ["position"], ascending = True, na_position = "last")
    df = df.reset_index(drop = True)
    df["Transcript Consequence"] = df["Transcript Consequence"].str.strip()
    return df
