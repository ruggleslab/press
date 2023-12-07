def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, label_binarize

from sklearn.decomposition import PCA, NMF
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import roc_auc_score, classification_report
from sklearn.metrics import roc_curve

from sklearn.base import clone
from joblib import load, dump
import random

sys.path.append('/Users/muller/Documents/RugglesLab')
from Pytools.plotting import plot_training_roc_curve_ci, plot_roc_curve
from Pytools.stats import bootstrap_auc_confidence, mean_confidence_interval

# All genes (not just PRESS)
path = 'data/clean/'
pace_all_genes = pd.read_csv('data/all_genes_normcounttab.txt', sep="\t", index_col=0).sort_index().T
pace_labels = pd.read_csv(path+'pace/labels.csv').to_numpy().ravel()

duke_all_genes = pd.read_csv('data/duke_validation_run3/normcounttab_genesymbols.csv', index_col=0).T.sort_index().T
duke_labels = pd.read_csv(path+'duke/labels_group1.csv').to_numpy().ravel()

aucs = []
for i in range(100):
    random_ints = random.sample(range(1, pace_all_genes.shape[1]), 402)
    # Model clone
    model = clone(load('models/jobs/rf-testing.joblib'))[-1]

    pace_random_genes = pace_all_genes.T.iloc[random_ints].T
    pace_random_genes = pace_random_genes.iloc[range(84)]


    duke_random_genes = duke_all_genes.iloc[:,random_ints]
    duke_random_genes = duke_random_genes.iloc[range(68)]

    model.fit(pace_random_genes, pace_labels)

    # print(duke_random_genes.head())
    print(f"Iteration Numer {i}")
    # ## Metrics!
    dump(model, 'models/jobs/null_model.joblib')

    y_pred = model.predict_proba(duke_random_genes)[:,1]
    score = roc_auc_score(duke_labels, y_pred)
    aucs.append(score)
    # bootstrap_auc_confidence(y_pred, duke_labels, n_bootstraps=1000)

working_roc_auc, confidence_lower, confidence_upper = mean_confidence_interval(aucs)
print(f"Confidence interval for the score: {working_roc_auc:0.2f} [{confidence_lower:0.2f} - {confidence_upper:0.2f}]")