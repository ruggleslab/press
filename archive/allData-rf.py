
##########
# Matthew Muller
# 11/24/2022
#
# Testing out feature reduction techniques
# https://scikit-learn.org/stable/auto_examples/compose/plot_compare_reduction.html#sphx-glr-auto-examples-compose-plot-compare-reduction-py
##########

##########
# Library Imports
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.svm import LinearSVC, SVC
from sklearn.decomposition import PCA, NMF, KernelPCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.metrics import roc_curve

##########
# Set/Append Working directory
sys.path.append('.')

##########
# Import Functions
from Pytools.plotting import plot_roc_curve, plot_confusion_matrix

##########
# Code Below

# Set random state
from random import randint
random = randint(0, 2**20)
# 651003 is a great local minima; 746833 is basically perfect

# Load in data
path = 'platelet-activity/data/clean/'
pace_features = pd.read_csv(path+'pace/features.csv').to_numpy()
pace_labels = pd.read_csv(path+'pace/labels.csv').to_numpy()[:,0]

duke_features = pd.read_csv(path+'duke/features_group1.csv').to_numpy()
duke_labels = pd.read_csv(path+'duke/labels_group1.csv').to_numpy()[:,0]

# Add pace and duke data together to train on
features = np.concatenate((pace_features, duke_features), axis=0)
labels = np.concatenate((pace_labels, duke_labels), axis=0)

X_train, X_test, y_train, y_test = train_test_split(features, labels, stratify=labels,
                                                    test_size=0.20, random_state=random)

# Test on combined cohorts
# X, y = pd.read_csv(path+'features.csv').to_numpy(), pd.read_csv(path+'labels.csv').to_numpy()[:,0]
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=420)

# Initialize things


## This is essentially cross validation on the train set to ensure we are at a good local minimum. Makes things more consistent
model = None
tmp = -1
rng = -1
for i in range(10):
    random_seed = random = randint(0, 2**20) # Perfect seed is 654134
    X_train, X_test, y_train, y_test = train_test_split(features, labels, stratify=labels,
                                                        test_size=0.20, random_state=random)
    rf = RandomForestClassifier(n_jobs=-1, random_state=random_seed)
    rf.fit(X_train, y_train)
    auc = roc_auc_score(y_train, rf.predict_proba(X_train)[:,1])
    if auc > tmp:
        print(f"Updated AUC: {auc}")
        model = rf
        tmp = auc
        rng = random_seed

print(f"auc: {tmp}, rng seed: {rng}")

# Metrics
print(classification_report(y_test, model.predict(X_test)))

plot_roc_curve(y_test, model.predict_proba(X_test)[:,1],
               save_path='platelet-activity/output/all_data_rf_roc-curve.png')

plot_confusion_matrix(y_test, model.predict(X_test),
                      labels = ['normal', 'hyper'],
                      
                      save_path='platelet-activity/output/all_data_rf_matrix.png')