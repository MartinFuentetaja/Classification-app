#Lo primero que debemos hacer es importar los modulos que necesitemos
##################################################################################################################################################################################
import numpy             as np
import matplotlib.pyplot as plt
import seaborn           as sns

import sklearn
from sklearn         import metrics
from sklearn.metrics import classification_report

##################################################################################################################################################################################

##################################################################################################################################################################################
#Esta función nos permite calcular la predicción con el modelo que le indiquemos.
def y_prediction(X_train, Y_train, x_test, y_test, model):
    model.fit(X_train, Y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test  = model.predict(x_test)
    score_train  = model.score(X_train, Y_train)
    score_test   = model.score(x_test, y_test)
    return y_pred_train, y_pred_test, score_train, score_test
#Nos devuelve la confusion matrix la cual nos permite analizar si el modelo es correcto dado que es posible ver una relación entre tanto con falsos positivos
#y negativos, como com los verdaderos.
def Confusion_Matrix(y_test, y_pred, score, path):
    conf_matrix = metrics.confusion_matrix(y_test, y_pred)
    path_mod = path + "_Confusion_Matrix"
    plt.figure(figsize=(8,8))
    sns.heatmap(conf_matrix, annot=True, fmt=".3f", linewidths=.5, square = True, cmap = 'Blues_r');
    plt.suptitle("Confusion Matrix")
    plt.ylabel('Actual label');
    plt.xlabel('Predicted label');
    all_sample_title = 'Accuracy Score: {0}'.format(score)
    plt.title(all_sample_title, size = 15);
    plt.savefig(path_mod)

#Esta nos permite hacer una gráfica de los datos realizando una clasificación entre los mismos.
def plot_decision_boundaries(X, y, model_class):
    """Function to plot the decision boundaries of a classification model.
    This uses just the first two columns of the data for fitting 
    the model as we need to find the predicted value for every point in 
    scatter plot.
    
    One possible improvement could be to use all columns fot fitting
    and using the first 2 columns and median of all other columns
    for predicting.
    
    Adopted from:
    """
    reduced_data = X[:, :2]
    model = model_class
    model.fit(reduced_data, y)

    # Step size of the mesh. Decrease to increase the quality of the VQ.
    h = .02     # point in the mesh [x_min, m_max]x[y_min, y_max].    

    # Plot the decision boundary. For that, we will assign a color to each
    x_min, x_max = reduced_data[:, 0].min() - 1, reduced_data[:, 0].max() + 1
    y_min, y_max = reduced_data[:, 1].min() - 1, reduced_data[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    # Obtain labels for each point in mesh using the model.
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()])    

    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.1),
                         np.arange(y_min, y_max, 0.1))

    Z = model.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)

    plt.contourf(xx, yy, Z, alpha=0.4)
    plt.scatter(X[:, 0], X[:, 1], c=y, alpha=0.8)
    plt.show()

def class_report_writting(y_test, y_pred, path):
    target_names = ["Benign", "Pathogenic"]
    file = open(path, "a")
    file.write(f"{classification_report(y_test, y_pred, target_names = target_names, digits = 6)} \n")
    file.close()
##################################################################################################################################################################################
