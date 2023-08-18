#Lo primero que debemos hacer es importar los modulos que necesitemos
##################################################################################################################################################################################

import pandas            as pd
import matplotlib.pyplot as plt
import numpy             as np
import seaborn           as sns
import os
import sys

import MachineLearningPackage.Directory_creation        as dc_create
import MachineLearningPackage.Tuning_Hyperparameters_CV as tuning_hyp_CV
import MachineLearningPackage.Sklearn_module            as sk_mod
import MachineLearningPackage.argparse_decision         as argparse

import sklearn 
from sklearn.preprocessing   import StandardScaler
from sklearn.linear_model    import LogisticRegression
from sklearn.neighbors       import KNeighborsClassifier
from sklearn.ensemble        import RandomForestClassifier
from sklearn.svm             import SVC
from sklearn.model_selection import train_test_split, cross_validate, cross_val_predict
from sklearn.metrics         import f1_score


##################################################################################################################################################################################

#Ahora debemos definir las funciones que vamos a utilizar. Cabe destacar que la mayoría están definidas en los módulos
#, usecols = ["My_Label", "A_i", "C_i", "D_i", "E_i", "F_i", "G_i", "H_i", "I_i"]
#["My_Label", "d_size", "d_vol", "d_pol_e", "d_ip_e", "d_hf_e", "d_msa", "d_hf","residue_conserv"]
#["My_Label", "d_size", "d_vol", "d_pol_e", "d_ip_e", "d_hf_e", "d_msa", "d_hf","residue_conserv", "pLDDT", "MyMTR_11", "MyMTR_15", "MyMTR_21", "MyMTR_31_C"]
#Lista = [["My_Label", "pLDDT", "MyMTR_11"], ["My_Label", "d_size", "d_vol"], ["My_Label", "d_size", "d_dpol_e"], ["My_Label", "d_size", "residue_conserv"],["My_Label", "d_size", "d_vol", "d_pol_e", "d_ip_e", "d_hf_e", "d_msa", "d_hf","residue_conserv"],["My_Label", "d_size", "d_vol", "d_pol_e", "d_ip_e", "d_hf_e", "d_msa", "d_hf","residue_conserv", "pLDDT", "MyMTR_11"], ["My_Label", "A_i", "C_i", "D_i", "E_i", "F_i", "G_i", "H_i", "I_i"]]
##################################################################################################################################################################################
def read_data_mod(file, columns):
    file_readed = pd.read_excel(file, header = 0, usecols = columns) #No hemos definido index.col para que se haga más fácil la lectura de los features. 
    #file_readed.drop(["Position", "Variant", "FlagRestrictivo", "FlagAmplio", "A_i", "C_i", "D_i", "E_i", "F_i", "G_i", "H_i", "I_i", "K_i", "L_i", "M_i", "N_i", "P_i", "Q_i", "R_i", "S_i", "T_i", "V_i", "W_i", "Y_i", "A_f", "C_f", "D_f", "E_f", "F_f", "G_f", "H_f", "I_f", "K_f", "L_f", "M_f", "N_f", "P_f", "Q_f", "R_f", "S_f", "T_f", "V_f", "W_f", "Y_f", "d_vol", "d_pol_e", "d_ip_e", "d_hf_e", "d_msa", "d_hf", "MyMTR_15", "MyMTR_21", "MyMTR_31_C"], axis = 1, inplace = True)
    file_readed = file_readed.iloc[:,:]  #:99 "MyMTR_11; :78 "na_to_na"
    return file_readed

#Obtener los datos reales para y. Estos serán fundamentales para poder predecirlos
def target_variables(dataFrame):
    y = np.array(dataFrame["My_Label"])
    if ((dataFrame[dataFrame["My_Label"] == 1].shape[0] == dataFrame.shape[0]) or (dataFrame[dataFrame["My_Label"] == 0].shape[0] == dataFrame.shape[0])): 
      sys.exit("It is not possible to do a Machine Learning classification for the indicated gene as there is a unique value to classify.")
    return y

#Obtener los Inputs o feature para poder crear el modelo. Solamente devuelve el feature x_i
def inputs_variables(dataFrame):
    data_mod = dataFrame.loc[:, dataFrame.columns != "My_Label"]
    rows, columns = data_mod.shape
    x = np.empty([rows, columns])
    for i in range(rows):
        x_i = np.array(data_mod.loc[i], ndmin = 2)
        x[i] = x_i
    return x

def sklearn_func(X_train_SK, Y_train, x_test_SK, y_test, model, path):
    path_test  = path + "_test_"
    path_train = path + "_train_"
    y_pred_SK_train, y_pred_SK_test, score_SK_train, score_SK_test = sk_mod.y_prediction(X_train_SK, Y_train, x_test_SK, y_test, model)
    sk_mod.Confusion_Matrix(Y_train, y_pred_SK_train, score_SK_train, path_train)
    sk_mod.Confusion_Matrix(y_test, y_pred_SK_test, score_SK_test, path_test)
    return y_pred_SK_train, y_pred_SK_test, score_SK_train, score_SK_test

def CrossValidation_analysis_RANDOM(model, grid_dictionary, number_splits, number_combinations, X_train, Y_train, x_test, y_test, path):
    Best_model_rand, random_score_rand, random_parameters_rand = tuning_hyp_CV.RandomSearch_Cross_Validation(model, grid_dictionary, number_splits, number_combinations, X_train, Y_train)
    y_pred_rand, score_rand = sklearn_func(X_train, Y_train, x_test, y_test, Best_model_rand, path)
    return Best_model_rand, random_score_rand, score_rand

def CrossValidation_analysis_GRID(model, grid_dictionary, number_splits, X_train, Y_train, x_test, y_test, path):
    Best_model_grid, score_grid_CV, parameters_grid_CV, Best_sandard_dev_train, Best_sandard_dev_test = tuning_hyp_CV.GridSearch_Cross_Validation(model, grid_dictionary, number_splits, X_train, Y_train)
    y_pred_grid_train, y_pred_grid_test, score_grid_train, score_grid_test = sklearn_func(X_train, Y_train, x_test, y_test, Best_model_grid, path)
    return Best_model_grid, score_grid_CV, parameters_grid_CV, Best_sandard_dev_train, Best_sandard_dev_test, y_pred_grid_train, y_pred_grid_test, score_grid_train, score_grid_test


def cross_validation(model, X_train, Y_train, n_splits):
    scoring             = ['accuracy', 'precision', 'recall', 'f1']
    results             = cross_validate(estimator = model, X = X_train, y = Y_train, cv = n_splits, scoring = scoring, verbose = 2, n_jobs = -1,  return_train_score = True)
    Train_Accuracy      = results['train_accuracy'].mean()
    Validation_Accuracy = results['test_accuracy'].mean()
    standard_deviation_train  = np.std(results['train_accuracy']).mean()
    standard_deviation_test   = np.std(results['test_accuracy']).mean()
    return Train_Accuracy, Validation_Accuracy, standard_deviation_train, standard_deviation_test              

#Vamos a definir una función para ejecutar los modelos y escribir los resultados

def execute_models_and_write_results(model_str, model, grid_CV, X_train, Y_train, x_test, y_test, number_splits, n_splits_CV, path, file, file_train, file_clss, score_train, score_test):
    path_grid = path + "_grid"
    path_CV   = path + "_CV_" + str(n_splits_CV)
    Best_model_grid, score_grid_CV, parameters_grid_CV, Best_sandard_dev_grid_train, Best_sandard_dev_grid_test, y_pred_grid_train, y_pred_grid_test, score_grid_train, score_grid_test = CrossValidation_analysis_GRID(model, grid_CV, number_splits, X_train, Y_train, x_test, y_test, path_grid)
    Train_Accuracy, Validation_Accuracy, standard_deviation_train, standard_deviation_test  = cross_validation(Best_model_grid, X_train, Y_train, n_splits_CV)
    y_pred_CV = cross_val_predict(estimator = Best_model_grid, X = X_train_mod, y = Y_train, cv = n_splits_CV)
    accuracy  = np.mean(y_pred_CV == Y_train)
    sk_mod.Confusion_Matrix(Y_train, y_pred_CV, accuracy, path_CV)
    #return parameters_grid_CV, score_train, score_grid_train, Best_sandard_dev_grid_train, score_test, score_grid_test, Best_sandard_dev_grid_test, n_splits_CV, Train_Accuracy, standard_deviation_train

    file.write(f"Results for {model_str}: \n")
    file.write(f"Best hyperparameters: {parameters_grid_CV} \n")
    file.write(f"Before Tuning train score: {float(score_train):8.5f} vs After tuning train score: {float(score_grid_train):8.5f} \u00B1 {float(Best_sandard_dev_grid_train):8.5f}. \n")
    file.write(f"Before Tuning test score: {float(score_test):8.5f} vs After tuning test score: {float(score_grid_test):8.5f} \u00B1 {float(Best_sandard_dev_grid_test):8.5f} \n")
    file.write(f"After Cross Validation (cv = {n_splits_CV}) train accuracy: {float(Train_Accuracy):8.5f} \u00B1 {float(standard_deviation_train):8.5} \n")
    file.write("\n")

    file_train.write(f"Results for {model_str}: \n")
    file_train.close()
    sk_mod.class_report_writting(Y_train, y_pred_CV, path_train)

    file_clss.write(f"Results for {model_str}: \n")
    file_clss.close()
    sk_mod.class_report_writting(y_test, y_pred_grid_test, path_clss)

##################################################################################################################################################################################

#############################################

#Cargamos los argumentos de argparse

args = argparse.parser_argument()

gen_name_upper = args.gen_name.upper()


general_path  = os.path.dirname(args.config_path)
Database_path = os.path.join(general_path, f"RESULTS/{args.gen_name}/{args.gen_name}_encoded.xlsx")

#Finalmente, debemos definir las variables

############################################
#Modelos:

standard_sc    = StandardScaler()
Logistic_Regr  = LogisticRegression()
KNearest       = KNeighborsClassifier()
randomForest   = RandomForestClassifier()
supportVector  = SVC()

##################################################################################################################################################################################

#Cargamos la información de la base de datos:
try:
   df             = read_data_mod(Database_path, args.column_names)
except FileNotFoundError:
   sys.exit()
x_train        = inputs_variables(df)
y_train        = target_variables(df)
rows, columns  = x_train.shape

#Creamos los paths para guardar la informacion
Results_path   = dc_create.control_directory(os.path.join(general_path, "RESULTS"))
gene_name_path = dc_create.control_directory(os.path.join(Results_path, args.gen_name))

directory, path = dc_create.creating_directory(df, gene_name_path)
graphs_path     = dc_create.control_directory(os.path.join(directory, "GRAPHS"))
graphs_path     = os.path.join(graphs_path, os.path.basename(path))

path_LR        = graphs_path + "_LR_"
path_KN        = graphs_path + "_KN_"
path_RF        = graphs_path + "_RF_"
path_SVC       = graphs_path + "_SVC_"
X_train, x_test, Y_train, y_test = train_test_split(x_train, y_train, test_size = 0.25, random_state = 0)


##############################################################################################################
#Aquí hacemos un feature scaling tanto del train_set como del test_set
X_train_mod = standard_sc.fit_transform(X_train)
x_test_mod  = standard_sc.transform(x_test) #Aquí no aplicamos fit_transform porque no queremos aplicar las misma transformaciones al train_set y al test_set
###################################################
#Lo primero que vamos a hacer es analizar las predicciones de cada modelo con los default hiperparametros y sin cross validation

#Con Logistic Regression con default hyperparameters:
y_pred_LR_train, y_pred_LR_test, score_LR_train, score_LR_test     = sklearn_func(X_train_mod, Y_train, x_test_mod, y_test, Logistic_Regr, path_LR)
#Con K-Nearest Neighbour con default hyperparameters:
y_pred_KN_train, y_pred_KN_test, score_KN_train, score_KN_test     = sklearn_func(X_train_mod, Y_train, x_test_mod, y_test, KNearest, path_KN)
#Con Random Forest con default hyperparameters:
y_pred_RF_train, y_pred_RF_test, score_RF_train, score_RF_test     = sklearn_func(X_train_mod, Y_train, x_test_mod, y_test, randomForest, path_RF)
#Con Support Vector Machine con default hyperparameters:
y_pred_SVC_train, y_pred_SVC_test, score_SVC_train, score_SVC_test = sklearn_func(X_train_mod, Y_train, x_test_mod, y_test, supportVector, path_SVC)
############################################################################################################
#Ahora vamos a construir los diccionarios con los hiperparametros más relevantes para cada caso

grid_LR_CV    = dict(C = list(np.arange(0.001, 10.0, 0.05)), solver = ["liblinear", "sag", "saga", "lbfgs"])
grid_KN_CV    = dict(metric = ["minkowski", "euclidean", "manhattan"], n_neighbors = list(range(1,20,1))) 
grid_RF_CV    = dict(max_depth = list(range(1,25,1)), max_features = ['sqrt','log2', None], n_estimators = list(range(1, 25, 1)))
grid_SVC_CV   = dict(kernel = ["linear", "poly", "rbf", "sigmoid"], C = list(np.arange(0.001, 200.0, 5.0)), gamma = list(np.arange(0.01,3.0,0.5)))

###########################################
#Ahora vamos a definir el número de splits que queremos para el grid Search

number_splits = int(args.number_splits)
n_splits_CV   = int(args.number_splits_grid)

###########################################

#Ahora vamos a definir un fichero donde guardaremos los datos

path_file  = path + "_" + str(n_splits_CV) + "_Best_models.txt"
path_clss  = path + "_test.txt"
path_train = path + "_" + str(n_splits_CV) + "_train.txt"
file       = open(path_file, "w")
file_clss  = open(path_clss, "w")
file_train = open(path_train,"w")

###########################################
#Con Logistic Regression:

#Una vez sabemos que C = 100.0 ahora debemos hacer un gridSearchCV para ver si podemos mejorar el score

execute_models_and_write_results("Logistic Regression", Logistic_Regr, grid_LR_CV, X_train_mod, Y_train, x_test_mod, y_test, number_splits, n_splits_CV, path_LR, file, file_train, file_clss, score_LR_train, score_LR_test)

"""
path_LR_grid = path_LR + "_grid"
path_LR_CV   = path_LR + "_CV_" + str(n_splits_CV)
Best_model_LR_grid, score_LR_grid_CV, parameters_LR_grid_CV, Best_sandard_dev_LR_grid_train, Best_sandard_dev_LR_grid_test, y_pred_LR_grid_train, y_pred_LR_grid_test, score_LR_grid_train, score_LR_grid_test = CrossValidation_analysis_GRID(Logistic_Regr, grid_LR_CV, number_splits, X_train_mod, Y_train, x_test_mod, y_test, path_LR_grid)
Train_Accuracy_LR, Validation_Accuracy_LR, standard_deviation_train_LR, standard_deviation_test_LR  = cross_validation(Best_model_LR_grid, X_train_mod, Y_train, n_splits_CV)
y_pred_LR_CV = cross_val_predict(estimator = Best_model_LR_grid, X = X_train_mod, y = Y_train, cv = n_splits_CV)
accuracy_LR  = np.mean(y_pred_LR_CV == Y_train)
sk_mod.Confusion_Matrix(Y_train, y_pred_LR_CV, accuracy_LR, path_LR_CV)
file.write("Results for Logistic Regression: \n")
file.write(f"Best hyperparameters: {parameters_LR_grid_CV} \n")
file.write(f"Before Tuning train score: {float(score_LR_train):8.5f} vs After tuning train score: {float(score_LR_grid_train):8.5f} \u00B1 {float(Best_sandard_dev_LR_grid_train):8.5f}. \n")
file.write(f"Before Tuning test score: {float(score_LR_test):8.5f} vs After tuning test score: {float(score_LR_grid_test):8.5f} \u00B1 {float(Best_sandard_dev_LR_grid_test):8.5f} \n")
file.write(f"After Cross Validation (cv = {n_splits_CV}) train accuracy: {float(Train_Accuracy_LR):8.5f} \u00B1 {float(standard_deviation_train_LR):8.5} \n")
file.write("\n")

file_train.write("Results for Logistic Regression: \n")
file_train.close()
sk_mod.class_report_writting(Y_train, y_pred_LR_CV, path_train)

file_clss.write("Results for Logistic Regression: \n")
file_clss.close()
sk_mod.class_report_writting(y_test, y_pred_LR_grid_test, path_clss)
"""
###########################################

###########################################
#Con K-Nearest:

#Una vez sabemos que n_neighbors = 11 ahora debemos hacer un gridSearchCV para ver si podemos mejorar el score
file_train = open(path_train, "a")
file_clss  = open(path_clss, "a")
execute_models_and_write_results("K-Nearest Neighbor", KNearest, grid_KN_CV, X_train_mod, Y_train, x_test_mod, y_test, number_splits, n_splits_CV, path_KN, file, file_train, file_clss, score_KN_train, score_KN_test)


"""
path_KN_grid = path_KN + "_grid"
path_KN_CV   = path_KN + "_CV_" + str(n_splits_CV)
Best_model_KN_grid, score_KN_grid_CV, parameters_KN_grid_CV, Best_sandard_dev_KN_grid_train, Best_sandard_dev_KN_grid_test, y_pred_KN_grid_train, y_pred_KN_grid_test, score_KN_grid_train, score_KN_grid_test = CrossValidation_analysis_GRID(KNearest, grid_KN_CV, number_splits, X_train_mod, Y_train, x_test_mod, y_test, path_KN_grid)
Train_Accuracy_KN, Validation_Accuracy_KN, standard_deviation_train_KN, standard_deviation_test_KN = cross_validation(Best_model_KN_grid, X_train_mod, Y_train, n_splits_CV)
y_pred_KN_CV = cross_val_predict(estimator = Best_model_KN_grid, X = X_train_mod, y = Y_train, cv = n_splits_CV)
accuracy_KN  = np.mean(y_pred_KN_CV == Y_train)
sk_mod.Confusion_Matrix(Y_train, y_pred_KN_CV, accuracy_KN, path_KN_CV)
file.write("Results for K-Nearest Neighbor: \n")
file.write(f"Best hyperparameters: {parameters_KN_grid_CV} \n")
file.write(f"Before Tuning train score: {float(score_KN_train):8.5f} vs After tuning train score: {float(score_KN_grid_train):8.5f} \u00B1 {float(Best_sandard_dev_KN_grid_train):8.5f}. \n")
file.write(f"Before Tuning test score: {float(score_KN_test):8.5f} vs After tuning test score: {float(score_KN_grid_test):8.5f} \u00B1 {float(Best_sandard_dev_KN_grid_test):8.5f}\n")
file.write(f"After Cross Validation (cv = {n_splits_CV}) train accuracy: {float(Train_Accuracy_KN):8.5f} \u00B1 {float(standard_deviation_train_KN):8.5f} \n")
file.write("\n")

file_train = open(path_train, "a")
file_train.write("Results for K-Nearest Neighbor: \n")
file_train.close()
sk_mod.class_report_writting(Y_train, y_pred_KN_CV, path_train)

file_clss = open(path_clss, "a")
file_clss.write("Results for K-Nearest Neighbor: \n")
file_clss.close()
sk_mod.class_report_writting(y_test, y_pred_KN_grid_test, path_clss)
"""
###########################################

###########################################
#Con Random Forest:

#RandomForestClassifier(max_depth=5, max_features='log2', n_estimators=10)

file_train = open(path_train, "a")
file_clss  = open(path_clss, "a")
execute_models_and_write_results("Random Forest", randomForest, grid_RF_CV, X_train_mod, Y_train, x_test_mod, y_test, number_splits, n_splits_CV, path_RF, file, file_train, file_clss, score_RF_train, score_RF_test)

"""
path_RF_grid = path_RF + "_grid"
path_RF_CV   = path_RF + "_CV_" + str(n_splits_CV)
Best_model_RF_grid, score_RF_grid_CV, parameters_RF_grid_CV, Best_sandard_dev_RF_grid_train, Best_sandard_dev_RF_grid_test, y_pred_RF_grid_train, y_pred_RF_grid_test, score_RF_grid_train, score_RF_grid_test = CrossValidation_analysis_GRID(randomForest, grid_RF_CV, number_splits, X_train_mod, Y_train, x_test_mod, y_test, path_RF_grid)
Train_Accuracy_RF, Validation_Accuracy_RF, standard_deviation_train_RF, standard_deviation_test_RF = cross_validation(Best_model_RF_grid, X_train_mod, Y_train, n_splits_CV)
y_pred_RF_CV = cross_val_predict(estimator = Best_model_RF_grid, X = X_train_mod, y = Y_train, cv = n_splits_CV)
accuracy_RF  = np.mean(y_pred_RF_CV == Y_train)
sk_mod.Confusion_Matrix(Y_train, y_pred_RF_CV, accuracy_RF, path_RF_CV)
file.write("Results for Random Forest: \n")
file.write(f"Best hyperparameters: {parameters_RF_grid_CV} \n")
file.write(f"Before Tuning train score: {float(score_RF_train):8.5f} vs After tuning train score: {float(score_RF_grid_train):8.5f} \u00B1 {float(Best_sandard_dev_RF_grid_train):8.5f}. \n")
file.write(f"Before Tuning test score: {float(score_RF_test):8.5f} vs After tuning test score: {float(score_RF_grid_test):8.5f} \u00B1 {float(Best_sandard_dev_RF_grid_test):8.5f} \n")
file.write(f"After Cross Validation (cv = {n_splits_CV}) train accuracy: {float(Train_Accuracy_RF):8.5f} \u00B1 {float(standard_deviation_train_RF):8.5f} \n")
file.write("\n")

file_train = open(path_train, "a")
file_train.write("Results for Random Forest: \n")
file_train.close()
sk_mod.class_report_writting(Y_train, y_pred_RF_CV, path_train)

file_clss = open(path_clss, "a")
file_clss.write("Results for Random Forest: \n")
file_clss.close()
sk_mod.class_report_writting(y_test, y_pred_RF_grid_test, path_clss)
"""

###########################################


###########################################
#Con Support Vector Machine:

#Ahora lo hacemos con grid

file_train = open(path_train, "a")
file_clss  = open(path_clss, "a")
execute_models_and_write_results("Support Vector Machine", supportVector, grid_SVC_CV, X_train_mod, Y_train, x_test_mod, y_test, number_splits, n_splits_CV, path_SVC, file, file_train, file_clss, score_SVC_train, score_SVC_test)

"""
path_SVC_grid = path_SVC + "_grid"
path_SVC_CV   = path_SVC + "_CV_" + str(n_splits_CV)
Best_model_SVC_grid, score_SVC_grid_CV, parameters_SVC_grid_CV, Best_sandard_dev_SVC_grid_train, Best_sandard_dev_SVC_grid_test, y_pred_SVC_grid_train, y_pred_SVC_grid_test, score_SVC_grid_train, score_SVC_grid_test = CrossValidation_analysis_GRID(supportVector, grid_SVC_CV, number_splits, X_train_mod, Y_train, x_test_mod, y_test, path_SVC_grid)
Train_Accuracy_SVC, Validation_Accuracy_SVC, standard_deviation_train_SVC, standard_deviation_test_SVC = cross_validation(Best_model_SVC_grid, X_train_mod, Y_train, n_splits_CV)
y_pred_SVC_CV = cross_val_predict(estimator = Best_model_SVC_grid, X = X_train_mod, y = Y_train, cv = n_splits_CV)
accuracy_SVC  = np.mean(y_pred_SVC_CV == Y_train)
sk_mod.Confusion_Matrix(Y_train, y_pred_SVC_CV, accuracy_SVC, path_SVC_CV)
file.write("Results for Support Vector Machine: \n")
file.write(f"Best hyperparameters: {parameters_SVC_grid_CV} \n")
file.write(f"Before Tuning train score: {float(score_SVC_train):8.5f} vs After tuning train score: {float(score_SVC_grid_train):8.5f} \u00B1 {float(Best_sandard_dev_SVC_grid_train):8.5f}. \n")
file.write(f"Before Tuning test score: {float(score_SVC_test):8.5f} vs After tuning test score: {float(score_SVC_grid_test):8.5f} \u00B1 {float(Best_sandard_dev_SVC_grid_test):8.5f} \n")
file.write(f"After Cross Validation (cv = {n_splits_CV}) train accuracy: {float(Train_Accuracy_SVC):8.5f} \u00B1 {float(standard_deviation_train_SVC):8.5f} \n")
file.close()

file_train = open(path_train, "a")
file_train.write("Results for Support Vector Machine: \n")
file_train.close()
sk_mod.class_report_writting(Y_train, y_pred_SVC_CV, path_train)

file_clss = open(path_clss, "a")
file_clss.write("Results for Support Vector Machine: \n")
file_clss.close()
sk_mod.class_report_writting(y_test, y_pred_SVC_grid_test, path_clss)
"""
##################################################

#print(f"Initial score for Logistic Regression {float(score_LR):8.5f} vs final score for Logistic Regression (RANDOM) {float(score_LR_rand):8.5f} vs final score for Logistic Regression (GRID) {float(score_LR_grid):8.5f}")
#print(f"Initial score for K-Nearest: {float(score_KN):8.5f} vs final score for K-Nearest (RANDOM): {float(score_KN_rand):8.5f} vs final score for K-Nearest (GRID): {float(score_KN_grid):8.5f}")


"""
Code_list                = ["My Code", "LR", "KN", "RF", "SVC"]
Accuracy_test_List_pLDDT = [95.714, 93.571, 95.714, 96.429, 93.571]
Accuracy_test_List_Structural = [94.286, 92.143, 94.286, 95.0, 95.714]
Accuracy_test_List_Integer = [91.429, 91.429, 85.714, 95.0, 93.571]
Accuracy_test_List_all = [94.286, 95.0, 93.571, 96.429, 94.286]

plt.grid()
plt.title("Models comparison")
plt.ylabel("Test accuracy(%)")
plt.xlabel("Models")
plt.ylim([80.0,100])
plt.bar(Code_list,Accuracy_test_List_Integer, width= 0.3, color = "r")
plt.show() 
"""
