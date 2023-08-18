#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd
import os
import requests
import time

#####################################################################################################################################################################################################################



def one_hot_amino(data):
    one_hot_encoded_data = pd.get_dummies(data, columns = ['initial amino acid', 'final amino acid'], dtype = int, prefix = ["Initial", "Final"], prefix_sep = "_")
    return one_hot_encoded_data

def one_hot_second_structure(data):
    one_hot_encoded_data = pd.get_dummies(data, columns = ['Second_Structure'], dtype = int, prefix = "", prefix_sep = "")
    headers = one_hot_encoded_data.columns.values.tolist()
    if (" " in headers):
        one_hot_encoded_data.rename(columns = {"H" : "Alpha helix", "B" : "Residue in isolated beta-bridge", "E" : "Extended strand",
                                               "G" : "3-helix", "I" : "5 helix", "T" : "hydrogen bonded turn", "S" : "bend",
                                               " " : "Loops and irregular elements"}, inplace = True)
    else:
        one_hot_encoded_data.rename(columns = {"C" : "coil", "H" : "helix", "E" : "beta_strand"}, inplace = True)
    return one_hot_encoded_data

def one_hot_general(data, column):
    one_hot_encoded_data = pd.get_dummies(data, columns = column, dtype = int, prefix = "", prefix_sep = "")
    return one_hot_encoded_data


def my_label_modification(df):
    conditions = (
        (df["Clinical significance (Last reviewed)"].str.split("(").str[0].str.lower().str.contains("benign|pathogenic"))
    )

    df_filtered = df[conditions].copy()
    df_filtered["My_Label"] = df_filtered["Clinical significance (Last reviewed)"].str.split("(").str[0].str.lower()
    df_filtered["My_Label"] = df_filtered["Clinical significance (Last reviewed)"].apply(lambda x: 0 if "benign" in x.split("(")[0].lower() else 1)
    
    my_label = df_filtered.pop("My_Label")
    df_filtered.insert(1, my_label.name, my_label)
    df_filtered.reset_index(drop = True)

    return df_filtered 


def data_Frame_modification(df):
    df["Name"] = df["Name"] + ", " + df["Protein change"]  + ", " + df["position"].astype(str)
    df.rename(columns = {'Name':'Variation'}, inplace = True)
    df.drop(columns = ["Gene(s)", "DNA change", "DNA position", "initial nucleotide", "final nucleotide", "Protein change (compl)","Protein change", "position", "Domain", "Function", "Region", "Condition(s)", "Review status", "Clinical significance (Last reviewed)", "Codons"], axis = 0, inplace = True)
    
    df = df.reset_index(drop = True)
    return df

