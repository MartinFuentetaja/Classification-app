
#Lo primero es importar los módulos que vamos a necesitar
##################################################################################################################################################################################
import sklearn
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
##################################################################################################################################################################################

#Aquí definimos las funciones que vamos a utilizar: RandomSearch_Cross_Validation, GridSearch_Cross_Validation. 
#Estas dos funciones nos permitirán hacer un análisis de la mejor combinación de hiperparametros para cierto módelo. De esta manera,
#seremos capaces de entrenar nuestro modelo con una mayor fiabilidad. 

#model               ---> Debemos introducir el modelo de clasificación que queremos analizar
#grid                ---> Un diccionario donde especificamos los hiperparametros a analizar
#number_splits       ---> Un integer donde indicamos como queremos que se divida el train_Set
#number_combinations ---> El número de combinaciones que queremos analizar
#X_train, Y_train    ---> Los train_Set con los que se entrenan al modelo

##################################################################################################################################################################################
def RandomSearch_Cross_Validation(model, grid, number_splits, number_combinations, X_train, Y_train):
    model_random = RandomizedSearchCV(estimator = model, param_distributions = grid, cv = number_splits, n_iter = number_combinations, n_jobs = -1, random_state = 0, scoring = "f1", error_score = "raise")
    model_random.fit(X_train, Y_train)
    random_score      = model_random.best_score_      #Devuelve el mejor score, es decir, para la mejor combinación el score que ha encontrado.
    Best_model        = model_random.best_estimator_  #Devuelve el modelo en principio con la mejor combinación de hiperparametros.
    random_parameters = model_random.best_params_     #Devuelve un diccionario con la mejor combinación de hiperparametros
    return Best_model, random_score, random_parameters

def GridSearch_Cross_Validation(model, grid, number_splits, X_train, Y_train, verbosity = 2):
    model_grid = GridSearchCV(estimator = model, param_grid = grid, cv = number_splits, n_jobs = -1, verbose = verbosity, scoring = "f1", error_score = "raise", return_train_score=True)
    model_grid.fit(X_train, Y_train)
    Best_score_CV      = model_grid.best_score_      #Devuelve la media cruzada del accuracy score de cada subset.
    Best_model         = model_grid.best_estimator_  #Devuelve el modelo en principio con la mejor combinación de hiperparametros.
    Best_parameters_CV = model_grid.best_params_     #Devuelve un diccionario con la mejor combinación de hiperparametros
    Best_sandard_dev_train  = model_grid.cv_results_['std_train_score'][model_grid.best_index_]
    Best_sandard_dev_test   = model_grid.cv_results_['std_test_score'][model_grid.best_index_]
    return Best_model, Best_score_CV, Best_parameters_CV, Best_sandard_dev_train, Best_sandard_dev_test
