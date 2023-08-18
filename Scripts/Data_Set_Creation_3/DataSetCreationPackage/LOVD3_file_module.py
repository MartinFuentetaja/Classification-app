#####################################################################################################################################################################################################################

import pandas as pd
import math
import re

#####################################################################################################################################################################################################################

"""
Esta primera función recoge el archivo Lovd3 en formato csv pero separado por tabulaciones. Su finalidad es realizar
una primera criva de los elementos de no interes: p.?, p.Ala872*...
Una vez hecho esto, devolverá un dataframe modificado sin dichos elementos.
"""

def modify_Lovd3(filename):
    df = pd.read_csv(filename, sep = "\t", header = 0, usecols = ["DNA change (cDNA)", "Protein", "Clinical classification"])
    df["Protein"] = df["Protein"].fillna("PASS")
    conditions = (
        (   
            ( 
                (~df["Protein"].str.contains("PASS")) &
                (~df["Protein"].str.contains(r"\=|\?|\*|\_|.0", na = False))
             ) &
            (
                (df["DNA change (cDNA)"].str.contains(">", na = False)) &
                (~df["DNA change (cDNA)"].str.contains("\+|\-|\_|\*", na = False))
             )
        ) &
            (df["Clinical classification"].str.lower().str.contains("likely|benign|pathogenic", na = False))
                  )
    df_filtered = df[conditions].copy()
    df_filtered["Protein"] = df_filtered["Protein"].replace(to_replace = "\(", value = "", regex = True).replace(to_replace = "\)", value = "", regex = True)
    df_filtered = df_filtered.reset_index(drop = True)
    return df_filtered

"""
Esta función recoge el dataframe anteriormente modificado y aplica una segunda criva para los elementos de la
columna "Clinical classification". En esta analizará los elementos repetidos, y de ser así los eliminará, es decir
los elementos en la columna de "Protein" que sean iguales y cuyos elementos de "Clinical classification" sean iguales.
"""

def second_modification_Lovd3(df):
    df["Duplicate"] = df.duplicated(subset = ["Protein", "DNA change (cDNA)"], keep = False)
    condition       = ( 
                    (df["Duplicate"]) &
                    (df["Clinical classification"].str.lower().str.contains("likely")
                        )
        )
    #Con esto de aquí dividimos los casos que contengan la palabra likely y que esten repetidos, y los que no.
    df_likely     = df[condition].drop_duplicates(subset = ["Protein", "DNA change (cDNA)"], keep = "first").reset_index(drop = True)
    df_not_likely = df[~condition].drop_duplicates(subset = ["Protein", "DNA change (cDNA)"], keep = "first").reset_index(drop = True)
    #Realizamos un .strip() por si hay algun carácter no visible que cause problemas al hacer el merge
    df_not_likely["DNA change (cDNA)"] = df_not_likely["DNA change (cDNA)"].str.strip()
    df_likely["DNA change (cDNA)"]     = df_likely["DNA change (cDNA)"].str.strip()
    df_likely["Protein"]               = df_likely["Protein"].str.strip()
    df_not_likely["Protein"]           = df_not_likely["Protein"].str.strip()
    #Hacemos un merge entre not_likely y likely para encontrar casos coincidentes, y de esta manera eliminarlos para poder unir los likely no repetidos con los not_likely.
    df_merge_in  = pd.merge(df_not_likely[["DNA change (cDNA)", "Protein", "Clinical classification"]], df_likely[["DNA change (cDNA)", "Protein", "Clinical classification"]], on = ["DNA change (cDNA)", "Protein"], how = "inner")
    df_merge_in  = df_merge_in.drop(columns = ["Clinical classification_x"]).rename(columns = {"Clinical classification_y" : "Clinical classification"})
    df_likely    = pd.merge(df_likely[["DNA change (cDNA)", "Protein", "Clinical classification"]], df_merge_in, on = ["DNA change (cDNA)", "Protein", "Clinical classification"], how = "outer", indicator = True)
    df_likely    = df_likely[df_likely['_merge'] == 'left_only'].reset_index(drop=True)
    
    df_merge     = pd.merge(df_not_likely[["DNA change (cDNA)", "Protein", "Clinical classification"]], df_likely[["DNA change (cDNA)", "Protein", "Clinical classification"]], on = ["DNA change (cDNA)", "Protein", "Clinical classification"], how = "outer")
    return df_merge

