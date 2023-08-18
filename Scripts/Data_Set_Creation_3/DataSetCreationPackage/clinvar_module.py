#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import pandas as pd
import re
import numpy as np
import sys

#####################################################################################################################################################################################################################

"""
Esta función recoge un archivo de texto descargado de clinvar. Hay que tener en cuenta que todos los archivos de clinvar tienen la misma estructura. Estan tabulados, y la primera linea corresponde al encabezado.
Lo primero leerá el archivo utilizando la librería pandas, y acto seguido creará una lista con los elementos de la columna "Name". De esta forma, comenzará a realizar una criva de los elementos que no le interesen. 
Sabiendo la estructura del fichero, los elementos de interes comienzan con "NM_". De esta manera, si el primer elemento no corresponde a "NM" lo saltará. Una vez haya encontrado una cadena que empiece por "NM" 
comenzará el segundo crivado dado que hay elementos que no son de interes.  Estos estarán compuestos en su mayoría por elementos de tipo: "*", "+", ">", "-". 
Así pues, una vez resuelta la criva, la función devolverá el DataFrame modificado. 
"""

def NM_modification_clinvar(txt_file):
    df = pd.read_csv(txt_file, sep = "\t", header = 0)
    conditions = ((df["Name"].str.contains("NM_")) & (~df["Name"].str.contains(r"\+|\-|\=|\*") & df["Name"].str.contains(">")))
    df_filtered = df[conditions].copy()
    df_filtered = df_filtered.reset_index(drop = True)
    if (df_filtered.empty):
       sys.exit("There are not missense variations for the indicated gene in ClinVar database.")
    return df_filtered

"""
Esta función recoge el DataFrame modificado anteriormente. También realiza modificaciones pero en la columna "Name" para poder extraer el cambio que se produce en la secuencia
del DNA. De esta forma, leerá cada línea de dicha columna y la información se la pasará a la función obtain_nucleotido_change(). Por lo tanto, esta devolverá tanto el cambio, 
la posición, y los nucleotidos inciales y finales, y guardará dicha información en respectivas listas. Una vez leídas todas las filas, se añadirán nuevas columnas para cada
nueva información, que se colocarán después de la columna "Gene(s)". Finalmente, la función devolverá el DataFrame modificado.
"""

def NM_second_modification_clinvar(df):
    rows, columns = df.shape
    nucl_change_list  = []
    nucl_pos_list     = []
    initial_nucl_list = []
    final_nucl_list   = []
    for i in range(rows):
        nucl_change, nucl_pos, initial_nucl, final_nucl = obtain_nucleotido_change(df.loc[i, "Name"])
        nucl_change_list.append(nucl_change)
        nucl_pos_list.append(nucl_pos)
        initial_nucl_list.append(initial_nucl)
        final_nucl_list.append(final_nucl)
    Name_loc = df.columns.get_loc("Gene(s)")
    df.insert(loc = Name_loc + 1, column = "DNA change", value = nucl_change_list)
    df.insert(loc = Name_loc + 2, column = "DNA position", value = nucl_pos_list)
    df.insert(loc = Name_loc + 3, column = "initial nucleotide", value = initial_nucl_list)
    df.insert(loc = Name_loc + 4, column = "final nucleotide", value = final_nucl_list)
    df_mod = df
    return df_mod

"""
Esta función recoge un DataFrame. Se dedica en hacer modificaciones en la columna "Protein change" utilizando la función extracting_protein_change().
De esta manera, leerá cada fila de la columna "Name", extraerá la posición, el aminoácido inicial,  aminoácido final, el cambio (letra_amino_incial + posición + letra_amino_final),
y el cámbio total, que se refiere al cambio pero indicando cada aminoácido con sus respectivas tres letras (e.j. Alanina -> "Ala"). De esta manera, comprueba si el aminoácido final
está con un "*", es decir, si corresponde a un "stop". De ser así, elimina dicha fila dado que no nos interesa. Finalmente, introduce toda la información en diferentes listas,
y las añade al DataFrame. Por lo tanto, el input es un DataFrame, y el output es el mismo DataFrame pero con nueva información de proteínas.
"""

def Protein_modification_clinvar_2(df_2):
    rows, columns     = df_2.shape
    position_list     = []
    init_amino_list   = []
    fin_amino_list    = []
    change_list       = []
    compl_change_list = []
    for i in range(rows):
        position, init_amino_mod, fin_amino_mod, change, compl_change = extracting_protein_change(df_2.loc[i, "Name"])
        if (fin_amino_mod == "*" or init_amino_mod == "*"):
            df_2 = df_2.drop(i)
        else:
            position_list.append(position)
            init_amino_list.append(init_amino_mod)
            fin_amino_list.append(fin_amino_mod)
            change_list.append(change)
            compl_change_list.append(compl_change)
    Protein_change_loc = df_2.columns.get_loc("Protein change")
    df_2 = df_2.drop(columns = ["Protein change"], axis = 1)
    df_2.insert(loc = Protein_change_loc, column = "Protein change (compl)", value = compl_change_list)
    df_2.insert(loc = Protein_change_loc + 1, column = "Protein change", value = change_list)
    df_2.insert(loc = Protein_change_loc + 2, column = "position", value = position_list)
    df_2.insert(loc = Protein_change_loc + 3, column = "initial amino acid", value = init_amino_list)
    df_2.insert(loc = Protein_change_loc + 4, column = "final amino acid", value = fin_amino_list)
    Review_status_loc  = df_2.columns.get_loc("Review status")
    df_2 = df_2.drop(df_2.iloc[:, Review_status_loc + 1:], axis = 1)
    df_2 = df_2.sort_values(by = ["position"], ascending = True, na_position = "last")
    df_2 = df_2.reset_index(drop = True)
    return df_2


"""
Esta función recoge una cadena de caracteres (e.j., NM_004519.4(KCNQ3):c.2537C>T (p.Thr846Met)), y su finalidad es extraer c.2537C>T, la cual refleja el cambio de nucleotidos en 
la cadena de DNA. De esta forma, el número inicial corresponde a la posición de la secuencia, "C" al nucleotido inical en esa posición y "T" al final. Esta información es fundamental
dado que nos permitirá calcular la variación que se de en la secuencia de la proteina, tanto en que posición como la variación de amino ácidos. Así pues, lo primero que hace es
un split en ":", por lo que devuelve "['NM_004519.4(KCNQ3)', 'c.2537C>T (p.Thr846Met)']".  Al solo interesarnos el elemnto situado en el índice [1], volveremos a hacer un split en "(".
Lo que ocurre es que al haber una tabulación, la suprimimos con sprit() antes de realizar el split("("). El output: ['c.2537C>T', 'p.Thr846Met)']. Hemos logrado nuestro objetivo
y lo guardamos en una lista. Lo que ocurre es que hacemos una subdivisión entre la posición y los nucleotidos. Para ello, realizamos un split en "." y acto seguido en ">", 
['2537C', 'T']. De aquí obtenemos, el nucleotido final. Por otro lado, del índice [0] extraemos la posición y el nucleotido inicial de la siguiente manera:  
nucl_pos     = int(re.findall('\d+', str_splt_4[0])[0]) y initial_nucl = str_splt_4[0][-1]. La función re.findall nos permite encontrar cualquier número en una cadena de texto, 
y nos devuelve una lista con los número encontrados. Finalmente el output es: ['c.2537C>T', 2537, 'C', 'T'], aunque no está en una lista.
"""

def obtain_nucleotido_change(string):
    str_splt     = string.split(":")
    str_splt_2   = str_splt[-1].strip().split("(")
    nucl_change  = str_splt_2[0]
    str_splt_3   = nucl_change.strip().split(".")
    str_splt_4   = str_splt_3[1].strip().split(">")
    nucl_pos     = int(re.findall('\d+', str_splt_4[0])[0])
    final_nucl   = str_splt_4[-1]
    initial_nucl = str_splt_4[0][-1]
    return nucl_change, nucl_pos, initial_nucl, final_nucl

"""
Esta función recoge una cadena de texto (e.j., NM_004519.4(KCNQ3):c.2537C>T (p.Thr846Met)) y devuelve "Thr846Met" desglosado como
position = 846, init_amino_mod = "Thr" y utilizando amino_Acid_classification() lo transformará en el amino ácido correspondiente. En este
caso "Thr" = "T". De la misma forma que devuelve el amino ácido final, y el cambio en la secuencia de la proteina "Thr846Met" a "T846M".
"""

def extracting_protein_change(string_chain):
    string_split_1  = string_chain.strip().split("(")
    string_change_1 = string_split_1[-1]
    compl_change    = string_split_1[-1][:-1]
    string_split_2  = string_change_1.strip().split(".")
    string_change_2 = string_split_2[-1]
    string_split_3  = string_change_2.strip().split(")")
    initial_amino   = string_split_3[0][:3]
    final_amino    = string_split_3[0][-3:]
    position       = int(re.findall('\d+', string_split_3[0])[0])
    init_amino_mod = amino_Acid_classification(initial_amino)
    fin_amino_mod  = amino_Acid_classification(final_amino)
    change         = init_amino_mod + str(position) + fin_amino_mod
    return position, init_amino_mod, fin_amino_mod, change, compl_change

"""
Esta función es muy simple. Recoge una cadena de texto cuya longitud es igual a 3, y corresponde a la clasificación de un aminoácido. De esta manera, devuelve la letra correspondiente
a dicho aminoácido. De no encontrarse entre los 20 aminoácidos comunes devuelve un "*".
"""

def amino_Acid_classification(amino_acid):
    if (amino_acid == "Ala"):
        amino_acid_mod = "A"
    elif (amino_acid == "Cys"):
        amino_acid_mod = "C"
    elif (amino_acid == "Asp"):
        amino_acid_mod = "D"
    elif (amino_acid == "Glu"):
        amino_acid_mod = "E"
    elif (amino_acid == "Phe"):
        amino_acid_mod = "F"
    elif (amino_acid == "Gly"):
        amino_acid_mod = "G"
    elif (amino_acid == "His"):
        amino_acid_mod = "H"
    elif (amino_acid == "Ile"):
        amino_acid_mod = "I"
    elif (amino_acid == "Lys"):
        amino_acid_mod = "K"
    elif (amino_acid == "Leu"):
        amino_acid_mod = "L"
    elif (amino_acid == "Met"):
        amino_acid_mod = "M"
    elif (amino_acid == "Asn"):
        amino_acid_mod = "N"
    elif (amino_acid == "Pro"):
        amino_acid_mod = "P"
    elif (amino_acid == "Gln"):
        amino_acid_mod = "Q"
    elif (amino_acid == "Arg"):
        amino_acid_mod = "R"
    elif (amino_acid == "Ser"):
        amino_acid_mod = "S"
    elif (amino_acid == "Thr"):
        amino_acid_mod = "T"
    elif (amino_acid == "Val"):
        amino_acid_mod = "V"
    elif (amino_acid == "Trp"):
        amino_acid_mod = "W"
    elif (amino_acid == "Tyr"):
        amino_acid_mod = "Y"
    else:
        amino_acid_mod = "*"
    return amino_acid_mod

