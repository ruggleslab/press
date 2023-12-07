
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import roc_auc_score
from sklearn.base import clone

from joblib import load

##########
# Set/Append Working directory
sys.path.append('/Users/muller/Documents/RugglesLab')

##########
# Import Functions
from MattTools.plotting import plot_training_roc_curve_ci
from MattTools.stats import bootstrap_auc_confidence

##########
# Code Below
path = 'data/clean/'
X_pace = pd.read_csv(path+'pace/features.csv')
y_pace = pd.read_csv(path+'pace/labels.csv').to_numpy()[:,0]

min_features_to_select = 1
# model = load('models/jobs/rf-testing.joblib')[-1] # testing model is a pipe!
model = load('models/jobs/extraTrees_trained_auc89.joblib')[-1], # this is a pipe, so I'm selecting just the model
cv = StratifiedKFold(5)

rfecv = RFECV(
    estimator=model,
    step=1,
    cv=cv,
    scoring="accuracy",
    min_features_to_select=min_features_to_select,
    n_jobs=-1,
)
rfecv.fit(X_pace, y_pace)

print(f"Optimal number of features: {rfecv.n_features_}")

# Graph our features
n_scores = len(rfecv.cv_results_["mean_test_score"])
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Mean test AUC")
plt.plot(
    range(min_features_to_select, n_scores + min_features_to_select),
    rfecv.cv_results_["mean_test_score"],
)
plt.fill_between(
        range(min_features_to_select, n_scores + min_features_to_select),
        rfecv.cv_results_["mean_test_score"]+rfecv.cv_results_["std_test_score"],
        rfecv.cv_results_["mean_test_score"]-rfecv.cv_results_["std_test_score"],
        color="grey",
        alpha=0.2,
        label="Standard Deviation",
    )
plt.axvline(rfecv.n_features_, linestyle='--', linewidth=1)
plt.title("Recursive Feature Elimination")
plt.show()