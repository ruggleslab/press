###########################################################################
#
#                            press_model
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: press_model

#======================== SETUP ========================#
#' The JSON object contains three key-value pairs:
#' * "normalization": This key corresponds to the method used for normalization. In this case, "mor" is used.
#' * "geneset": This key corresponds to the path of the selected features file. In this case, the path is "output/feature_selection/rfe/logreg/selected_features.csv".
#' * "outdir": This key corresponds to the output directory where the results will be stored. In this case, the output directory is "output/reduce_press/rfeLOGREG_mor/".
import argparse
args = argparse.ArgumentParser()
args.add_argument('--json', type=str, help='Path to the JSON input file')
args = args.parse_args()

import json
with open(args.json, 'r') as f:
    params = json.load(f)

# General libaries
import os
import sys
import time
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import dump, load
import plotnine as p9
# Modeling libraries
from sklearn import base
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import metrics
from sklearn import feature_selection

# Custom libraries
from MattTools import utils, plotting
## Set up directories
# Root directory
root_dir = os.getcwd()
# Output directory
experiment = "modeling"
indir = os.path.join(root_dir, params['outdir'], 'datasets') + "/"
outdir = os.path.join(root_dir, params['outdir'], experiment) + "/"
print(params['outdir'])
# Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

## Set up some general items
utils.set_random_seed(420)
utils.hide_warnings()

def get_data(indir, dataset):
    # check if the counts and labels exist
    if not os.path.exists(f'{indir}/{dataset}/normalized_counts.csv'):
        raise ValueError(f'No normalized counts for {dataset}')
    if not os.path.exists(f'{indir}/{dataset}/label.csv'):
        raise ValueError(f'No label for {dataset}')
    
    X = pd.read_csv(f'{indir}/{dataset}/normalized_counts.csv', index_col=0, header=0).T
    y = pd.read_csv(f'{indir}/{dataset}/label.csv', index_col=0, header=0).values
    return X, y

def get_model(json):
    # get the model from the json
    model = json['model']
    if model == 'PRESS451':
        m = VotingClassifier(
            estimators=
            [('rf'+str(i), RandomForestClassifier(class_weight='balanced', n_estimators=100)) for i in range(n)] +
            [('extraTrees'+str(i), ExtraTreesClassifier(class_weight='balanced', random_state=seeds[i])) for i in range(n) ] +
            [('gbc'+str(i), GradientBoostingClassifier(random_state=seeds[i], max_features='log2', n_estimators=60)) for i in range(n) ] +
            [('ada'+str(i), AdaBoostClassifier(learning_rate=10)) for i in range(n) ],
            voting='soft',n_jobs=-1,
        )
        parameters = {
            'rf0__n_estimators': [20, 50, 100],
            'rf0__class_weight': ['balanced', 'balanced_subsample'],
            'rf0__random_state': seeds,
            'extraTrees0__n_estimators': [10, 100, 200],
        }
    elif model == 'RF':
        m = RandomForestClassifier()
        parameters = {
            'n_estimators': [10, 100, 200, 500],
            'class_weight': ['balanced', 'balanced_subsample'],
        }
    elif model == 'ADA':
        m = AdaBoostClassifier()
        parameters = {
            'n_estimators': [10, 100, 200, 500],
            'learning_rate': [0.01, 0.1, 1, 10],
        }
    elif model == 'LOGREG':
        m = LogisticRegression()
        parameters = {
            'C': [0.01, 0.1, 1, 10, 100],
            'class_weight': ['balanced', None],
        }
    elif model == 'RUS':
        m = RUSBoostClassifier()
        parameters = {
            'n_estimators': [10, 50, 100, 200, 500],
            'learning_rate': [0.01, 0.1, 1, 10],
        }
    else:
        raise ValueError(f'Unknown model: {model}')
    return m, parameters

def test_model(X, y, outdir):
    os.makedirs(outdir, exist_ok=True)
    X = X[genes]
    X = pd.DataFrame(scaler.transform(X), index=X.index, columns=X.columns)
    
    report = metrics.classification_report(y, model.predict(X))
    with open(outdir+'classification_report.txt', 'w') as f:
        f.write(report)
        
    plotting.plot_roc_curve_ci(model, X, y, title='ROC Curve', save_path=outdir+'roc_curve.png')
    plotting.plot_prc_curve(model, X, y, title='PRC Curve', save_path=outdir+'prc_curve.png')
    plotting.plot_confusion_matrix(y, model.predict(X), outdir+'confusion_matrix.png')
    plt.close('all')
    
    # save the predictions, prediction probabilities, and true labels
    df = pd.DataFrame({
        'True': y.flatten(),
        'Predicted': model.predict(X),
        'Probabilities': model.predict_proba(X)[:,1]
    })
    sns.histplot(df['Probabilities'], kde=True)
    plt.savefig(outdir+'probabilities_hist.png')
    plt.close('all')
    sns.boxplot(x='True', y='Probabilities', data=df)
    plt.savefig(outdir+'probabilities_boxplot.png')
    df.to_csv(outdir+'predictions.csv', index=False)
    
    # make a summary of the groups
    summary = pd.DataFrame({
        'Dataset': [outdir.split('/')[-2]],
        'AUC': [metrics.roc_auc_score(y, model.predict_proba(X)[:,1])],
        'Accuracy': [metrics.accuracy_score(y, model.predict(X))],
        'Balanced Accuracy': [metrics.balanced_accuracy_score(y, model.predict(X))],
        'F1': [metrics.f1_score(y, model.predict(X))],
        'Precision': [metrics.precision_score(y, model.predict(X))],
        'Recall': [metrics.recall_score(y, model.predict(X))],
    })
    summary.to_csv(outdir+'summary.csv')
    return summary
    

#======================== CODE ========================#
# start off by running the datasets.R file
# os.system('Rscript code/01_prep_datasets.R' + ' ' + args.json)

genes = pd.read_table(params['geneset']).values.flatten()

genes_df = pd.DataFrame(genes, columns=['Gene'])
genes_df.to_csv(outdir + 'genes.csv', index=False)

X_pace, y_pace = get_data(indir, 'derivation')
X_pace = X_pace[genes]

scaler = preprocessing.StandardScaler()
X_pace = pd.DataFrame(scaler.fit_transform(X_pace), index=X_pace.index, columns=X_pace.columns)
dump(scaler, outdir+'scaler.joblib')

# print the shapes
print("X_pace shape:", X_pace.shape)
print("y_pace shape:", y_pace.shape)

#======================== Modeling ========================#
os.makedirs(outdir+'/modeling', exist_ok=True)
#### Make model pipeline (if needed) and search for params
# The general idea here is fit each gene to a SVC and then add them into a voting classifier
from sklearn.ensemble import VotingClassifier, StackingClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression, Perceptron
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from imblearn.ensemble import BalancedRandomForestClassifier, BalancedBaggingClassifier, RUSBoostClassifier

# good RF random seeds:
seeds = [1148093, 1095286, 1665788, 97057, 152878, 4543, 277452, 295106, 191278, 701043, 81388, 951209, 1327001, 527903, 1148093, 1095286, 1665788, 97057, 152878, 4543, 277452, 295106, 191278, 701043, 81388, 951209, 1327001, 527903]
n = 28
# parameters = {} # this takes too long let's just use the defaults

# make a voting classifier made a bunch of RFs
model, parameters = get_model(params)
cv = model_selection.GridSearchCV(
    model, parameters,
    n_jobs=-1,
    scoring="balanced_accuracy",
    )
cv.fit(X_pace, y_pace)
cv_results = cv.cv_results_
np.save(outdir+'/modeling/'+'cv_results', cv_results)

model = cv.best_estimator_
model.fit(X_pace, y_pace)
classification_report = metrics.classification_report(y_pace, model.predict(X_pace))
with open(outdir+'/modeling/'+'classification_report.txt', 'w') as f:
    f.write(classification_report)
dump(model, outdir+'model.joblib')
cv = model_selection.StratifiedKFold(n_splits=10)
plotting.plot_training_roc_curve_ci(
    model, X_pace, y_pace,
    title='10-fold Cross Validation Training on PACE Cohort',
    cv = cv,
    save_path = outdir+'/modeling/'+'pace_training_roc_curve.png'
    )


#======================== Validation ========================#
######## Duke group 1
os.makedirs(outdir+'/evaluation/', exist_ok=True)
# get the normalized counts and labels from the subdirectories
# make an empty report
report_all = pd.DataFrame()
for subdir in os.walk(indir).__next__()[1]:
    X, y = get_data(indir, subdir)
    report = test_model(X, y, outdir+'/evaluation/'+subdir+'/')
    report_all = report_all.append(report)
# write the dict of dicts out
with open(outdir+'/evaluation/'+'classification_reports.json', 'w') as f:
    json.dump(report_all, f, indent=4)