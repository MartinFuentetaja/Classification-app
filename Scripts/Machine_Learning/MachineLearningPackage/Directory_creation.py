##################################################################################################################################################################################

import os

##################################################################################################################################################################################

##################################################################################################################################################################################

#Esta función nos creará un directorio (en caso de que no exista) y almacernará todos los datos en el mismo. 
def creating_directory(df, Results_path):
    file_columns = df.columns.tolist()[1:] #Imponiendo que sea desde [1:] evitamos la columna que corresponde a "MyLabel"
    if (len(file_columns) <= 3):
        directory_name = "_vs_".join(file_columns) + "_data"
        file_name      = "_vs_".join(file_columns)
    else:
        directory_name = "_vs_".join(file_columns[0:2]) + "_data_" + str(len(file_columns)) + "_columns"
        file_name      = "_vs_".join(file_columns[0:2])
    
    #.join permite unir cadenas de texto para formar un path de manera correcta. \cwd\directory_name\
    directory_path = control_directory(os.path.join(Results_path, directory_name))
    path           = os.path.join(directory_path, file_name)
    return directory_path, path

#Controlar si existe un path:

def control_directory(directory_path):
    if not os.path.exists(directory_path):
        os.mkdir(directory_path)
    return directory_path

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
##################################################################################################################################################################################
