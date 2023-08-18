#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import numpy as np
import pandas as pd 

#####################################################################################################################################################################################################################

"""
Esta primera función recoge el dataFrame procedente de la unión entre df_codonProtein, df_GenomAD y df_CodigoGenetico. De esta manera,
tenemos un dataFrame que relaciona Codones con su respectivos amino ácidos por posición, y a su vez, la variación tanto en la proteína como en 
el DNA correspondiente a dicha posición. Además, contiene tanto el Allele Count, Number, Missense y Synonimous por posición. De esta manera,
esta función calculará la frecuencia sumando los Allele Count y Number del V2 y V3 entresí, y los dividirá, para volver a dividirlos por el número
de MissenseNotPoss en la próxima función. Una vez realizado esto, añade la columna de la frecuencia al dataFrame.
"""

def Frequency_calculation(df, length):
    df            = df.fillna({"Allele Count_x" : 0, "Allele Count_y" : 0, "Allele Number_x" : 0, "Allele Number_y" : 0})
    allele_count  = np.array(df[["Allele Count_x", "Allele Count_y"]].sum(axis = 1))
    allele_number = np.array(df[["Allele Number_x", "Allele Number_y"]].sum(axis = 1))
    frequency     = np.where(allele_number == 0, 0.0, np.nan_to_num(np.divide(allele_count, allele_number), nan = 0.0, posinf = 0.0)) #/ np.array(df["SynonymousPos"])
    df.insert(loc = df.columns.get_loc("SynonymousPos") + 1, column = "Frequency", value = frequency)
    MTI           = MTR_calculation(df, length, const = 1.0e-7)
    return MTI

"""
Esta función recoge el dataFrame resultante de la anterior función, y modifica el cáclculo para extraer un MTR correspondiente a la longitud de
la secuencia, es decir, a pesar de que por cada posición podamos llegar a tener diferentes variaciones, nos interesa tener un único valor por
posición. Por lo tanto, debemos recoger todos las frecuancias calculadas para la misma posición, sumarlas, y dividirlas por el valor correspondiente
de MissNoStopPos. De esta manera, obtendremos un único valor por posición dado que dividimos las frecuencias obtenidas por todas las posibles
frecuencias que se podrían haber dado en dicha posición. Este último valor es precisamente el que indica MissNoStopPos.
"""

def MTR_calculation(df, length, const = 1.0e-7):
    MTR = np.array([])
    for i in range(length):
        df_mod = df[["MissNoStopPos", "Frequency"]][df["position"] == i + 1].reset_index(drop=True)
        if (df_mod.empty):
            MTR = np.append(MTR, np.array([0]))
        else:
            number = df_mod["Frequency"].sum(axis = 0) / df_mod.loc[0, "MissNoStopPos"].astype(int)
            MTR = np.append(MTR, np.array([number]))
    MTR = MTR / (MTR + const)
    return MTR

"""
Esta función recoge el MTR obtenido de la anterior función, además de "promedio_number". Su finalidad es calcular el MTI, es decir,
el MTR ponderado. Lo que hace es relacionar los valores del MTR por posición dependiendo el valor introducido para "promedio_number", hasta
completar 31 ciclos. De esta manera, si promedio_number = 3, el MTI será el resultado de haber relacionado tres posiciones del MTR entresí. 
Al hacer esto, en el MTI tendremos por cada posición una media de diferentes posiciones del MTR.
"""

def MTI_ponderado(MTR, promedio_number):
    number  = int(promedio_number - 1)
    for i in range(1, 33, number):
        MTR_mod = np.array([])
        for j in range(MTR.shape[0]):
            if (j - int(number/2) < 0):
                MTR_mod = np.append(MTR_mod, np.mean(MTR[:number + 1]))
            elif (j == MTR.shape[0]):
                MTR_mod = np.append(MTR_mod, np.mean(MTR[j - int(number/2):]))
            else:
                MTR_mod = np.append(MTR_mod, np.mean(MTR[j - int(number/2):int(number/2) + j + 1]))
        MTR = MTR_mod
    return MTR_mod


