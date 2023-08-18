#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd


#####################################################################################################################################################################################################################



def polar_change(init_polar, fin_polar):
    if (init_polar == "non_polar" and fin_polar == "non_polar"):
        change = "np_to_np"
    elif (init_polar == "non_polar" and fin_polar == "polar"):
        change = "np_to_p"
    elif (init_polar == "polar" and fin_polar == "non_polar"):
        change = "p_to_np"
    elif (init_polar == "polar" and fin_polar == "polar"):
        change = "p_to_p"
    return change

def aromatic_change(init_arom, fin_arom):
    if (init_arom == "non_aromatic" and fin_arom == "non_aromatic"):
        change = "na_to_na"
    elif (init_arom == "non_aromatic" and fin_arom == "aromatic"):
        change = "na_to_a"
    elif (init_arom == "aromatic" and fin_arom == "non_aromatic"):
        change = "a_to_na"
    elif (init_arom == "aromatic" and fin_arom == "aromatic"):
        change = "a_to_a"
    return change

"""
Esta función recoge el amino ácido inical y final del dataFrame que estamos construyendo, y el dataFrame (que es fijo) de la lista de Amino ácidos.
De esta manera, hace una comparación entre el amino ácido final e incial de diferentes propiedades:

Categóricos:

charge_change -> Cambio en la carga (Si pasa de neutro a positivo: neu_to_pos. O de negativo a positivo: neg_to_pos...)
polarity_change -> Cambio en la polaridad (No polar a polar : np_to_p. Polar a polar: p_to_p ...)
aroma_change -> Cambio en la aromaticidad (No aromático a aromático: na_to_a ...)

Numéricas:

d_size -> Cambio en la carga
d_vol -> Cambio en el volumen
d_pol_e -> Cambio en la polarizabilidad
d_ip_e -> Cambio en el punto isométrico
d_hf_e -> Cambio en la hidrofobidad
d_msa - > Cambio en la "Mean Solvent Accesible Surface Area"

"""

def amino_acid_comparision(initial_amino_acid, final_amino_acid, df_amino):
    amino_Acid_features = df_amino
    init_amino       = amino_Acid_features.loc[initial_amino_acid]
    fin_amino        = amino_Acid_features.loc[final_amino_acid]
    charge_change    = init_amino["Charge"][:3] + "_to_" + fin_amino["Charge"][:3]
    polarity_change  = polar_change(init_amino["Polarity"], fin_amino["Polarity"])
    aroma_change     = aromatic_change(init_amino["Aromaticity"], fin_amino["Aromaticity"])
    d_size           = fin_amino["Molecular Weight"] - init_amino["Molecular Weight"]
    d_vol            = fin_amino["Volume"] - init_amino["Volume"]
    d_pol_e          = fin_amino["Polarizability"] - init_amino["Polarizability"]
    d_ip_e           = fin_amino["Isometric Point"] - init_amino["Isometric Point"]
    d_hf_e           = fin_amino["Hidrophobicity"] - init_amino["Hidrophobicity"]
    d_msa            = fin_amino["Mean Solvent Accesible Surface Area"] - init_amino["Mean Solvent Accesible Surface Area"]
    return charge_change, polarity_change, aroma_change, d_size, d_vol, d_pol_e, d_ip_e, d_hf_e, d_msa
