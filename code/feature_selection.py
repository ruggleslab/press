###########################################################################
#
#                            feature_selection
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-24
# Script Name: feature_selection

'''
The purpose of this script is to perform pythonic feature selection on the data. This will be done using a variety of methods, including Recursive Feature Elimination (RFE) and Sequential Feature Selection (SFS). The goal is to reduce the number of features to a manageable level, while still retaining the most important features. This will be done using a variety of classifiers, including Logistic Regression, Random Forest, and AdaBoost. The results will be saved to the output directory for further analysis.
'''


#======================== SETUP ========================#

# General libaries
import os
import sys
import time
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import load, dump
# Modeling libraries
from sklearn import base
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import metrics
# Custom libraries
from MattTools import utils
sys.path.append('.')
import src.data_functions as dfs
import src.feature_selection as fs
## Set up directories
# Root directory
root_dir = os.getcwd()
# Output directory
experiment = "feature_selection"
outdir = os.path.join(root_dir, "output", experiment) + "/"
# Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

## Set up some general items
utils.set_random_seed(420)
utils.hide_warnings()

#======================== Data Prep ========================#
X, y = dfs.load_data('data/clean/pace/features.csv', 'data/clean/pace/labels.csv')

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)
dump(scaler, 'models/press_reduction_scaler.joblib')
print('Original shape:', X.shape)

#======================== RFE ========================#
# So let's try out a few RFE classifiers and see how these are able to reduce the number of features
print('Running RFE')
outdir_rfe = outdir + "rfe/"
os.makedirs(outdir_rfe, exist_ok=True)

from sklearn.linear_model import LogisticRegression
print('Using Logistic Regression')
rfe_selector = fs.RFESelector(LogisticRegression())
rfe_selector.fit(X, y)
rfe_selector.plot_features(outdir_rfe+'logreg')
rfe_selector.save_results(outdir_rfe+'logreg')
dump(rfe_selector, outdir_rfe+'logreg' + 'rfe.joblib')

from sklearn.ensemble import RandomForestClassifier
print('Using Random Forest')
rfe_selector = fs.RFESelector(RandomForestClassifier())
rfe_selector.fit(X, y)
rfe_selector.plot_features(outdir_rfe+'rf')
rfe_selector.save_results(outdir_rfe+'rf')
dump(rfe_selector, outdir_rfe+'rf' + 'rfe.joblib')

from sklearn.ensemble import AdaBoostClassifier
print('Using AdaBoost')
rfe_selector = fs.RFESelector(AdaBoostClassifier())
rfe_selector.fit(X, y)
rfe_selector.plot_features(outdir_rfe+'ada')
rfe_selector.save_results(outdir_rfe+'ada')
dump(rfe_selector, outdir_rfe+'ada' + 'rfe.joblib')

#======================== SFS ========================#
# Sequential Feature Selection (SFS) is another method that can be used to select features. This also let's us try out PRESS451 again and see how it performs for minimal selection. Quite taxing O-wise though.
print('Running SFS')
outdir_sfs = outdir + "sfs/"
os.makedirs(outdir_sfs, exist_ok=True)

press451 = load('models/jobs/press451.joblib')
print('Using PRESS451')
SeqSelector = fs.SFSSelector(press451, 0.05)
SeqSelector.fit(X, y)
SeqSelector.plot_results(outdir_sfs+'press451')
SeqSelector.save_results(outdir_sfs+'press451')
dump(SeqSelector, outdir_sfs+'press451' + 'sfs.joblib')