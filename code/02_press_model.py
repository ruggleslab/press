###########################################################################
#
#                            press_model
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: press_model

#======================== SETUP ========================#
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

def test_model(X, y, outdir):
    os.makedirs(outdir, exist_ok=True)
    X = X[genes]
    X = pd.DataFrame(scaler.transform(X), index=X.index, columns=X.columns)
    
    report = metrics.classification_report(y, model.predict(X))
    with open(outdir+'classification_report.txt', 'w') as f:
        f.write(report)
        
    plotting.plot_roc_curve_ci(model, X, y, title='ROC Curve', save_path=outdir+'roc_curve.png')
    plotting.plot_confusion_matrix(y, model.predict(X), outdir+'confusion_matrix.png')
    plt.close('all')
    
    # save the predictions, prediction probabilities, and true labels
    df = pd.DataFrame({
        'True': y,
        'Predicted': model.predict(X),
        'Probabilities': model.predict_proba(X)[:,1]
    })
    sns.histplot(df['Probabilities'], kde=True)
    plt.savefig(outdir+'probabilities_hist.png')
    plt.close('all')
    sns.boxplot(x='True', y='Probabilities', data=df)
    plt.savefig(outdir+'probabilities_boxplot.png')
    df.to_csv(outdir+'predictions.csv', index=False)

#======================== CODE ========================#
# start off by running the datasets.R file
# os.system('Rscript code/01_prep_datasets.R' + ' ' + args.json)

genes = pd.read_table(params['geneset']).values.flatten()

genes_df = pd.DataFrame(genes, columns=['Gene'])
genes_df.to_csv(outdir + 'genes.csv', index=False)

X_pace = pd.read_csv(f'{indir}/derivation/normalized_counts.csv', index_col=0, header=0).T
X_pace = X_pace[genes]
y_pace = pd.read_csv(f'{indir}/derivation/label.csv')['hypercohort_inrnaseq_AP']

scaler = preprocessing.StandardScaler()
X_pace = pd.DataFrame(scaler.fit_transform(X_pace), index=X_pace.index, columns=X_pace.columns)
dump(scaler, outdir+'scaler.joblib')

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

# good RF random seeds:
seeds = [1148093, 1095286, 1665788, 97057, 152878, 4543, 277452, 295106, 191278, 701043, 81388, 951209, 1327001, 527903, 1148093, 1095286, 1665788, 97057, 152878, 4543, 277452, 295106, 191278, 701043, 81388, 951209, 1327001, 527903]
n = 28
parameters = {
    # 'rf0__max_features': ['log2', 'sqrt', 'auto'],
    'rf0__n_estimators': [20, 50, 100],
    'rf0__class_weight': ['balanced', 'balanced_subsample'],
    'rf0__random_state': seeds,
    'extraTrees0__n_estimators': [10, 100, 200],
}
parameters = {} # this takes too long let's just use the defaults

# make a voting classifier made a bunch of RFs
model = VotingClassifier(
    estimators=
    [('rf'+str(i), RandomForestClassifier(class_weight='balanced', n_estimators=100)) for i in range(n)] +
    
    [('extraTrees'+str(i), ExtraTreesClassifier(class_weight='balanced', random_state=seeds[i])) for i in range(n) ] +
    
    [('gbc'+str(i), GradientBoostingClassifier(random_state=seeds[i], max_features='log2', n_estimators=60)) for i in range(n) ] +
    
    [('ada'+str(i), AdaBoostClassifier(learning_rate=10)) for i in range(n) ],
    
    voting='soft',
    n_jobs=-1,
    )
# model = RandomForestClassifier(class_weight='balanced', n_estimators=100, random_state=seeds[0])

cv = model_selection.GridSearchCV(
    model, parameters,
    n_jobs=-1,
    scoring="roc_auc",
    )
cv.fit(X_pace, y_pace)
cv_results = cv.cv_results_
np.save(outdir+'/modeling/'+'cv_results.csv', cv_results)

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
X_duke_t1 = pd.read_csv(f'{indir}/duke_t1/normalized_counts.csv', index_col=0, header=0).T
y_duke_t1 = pd.read_csv(f'{indir}/duke_t1/label.csv')['hypercohort']
test_model(X_duke_t1, y_duke_t1, outdir+'/evaluation/'+'duke_t1/')

######## Duke group 2
X_duke_t2 = pd.read_csv(f'{indir}/duke_t2/normalized_counts.csv', index_col=0, header=0).T
y_duke_t2 = pd.read_csv(f'{indir}/duke_t2/label.csv')['hypercohort']
test_model(X_duke_t2, y_duke_t2, outdir+'/evaluation/'+'duke_t2/')

########## HARP MI_v_Obstr
X_harp = pd.read_csv(f'{indir}/harp/normalized_counts.csv', index_col=0, header=0).T
y_harp = pd.read_csv(f'{indir}/harp/label.csv')['MI_v_Obstr']
test_model(X_harp, y_harp, outdir+'/evaluation/'+'harp/')

########## SLE
X_sle = pd.read_csv(f'{indir}/sle/normalized_counts.csv', index_col=0, header=0).T
y_sle = pd.read_csv(f'{indir}/sle/label.csv')['Diagnosis']
test_model(X_sle, y_sle, outdir+'/evaluation/'+'sle/')

########## COVID
X_covid = pd.read_csv(f'{indir}/covid/normalized_counts.csv', index_col=0, header=0).T
y_covid = pd.read_csv(f'{indir}/covid/label.csv')['Cohort']
test_model(X_covid, y_covid, outdir+'/evaluation/'+'covid/')

########## MACLE2
X_pace = pd.read_csv(f'{indir}/MACLE2/normalized_counts.csv', index_col=0, header=0).T
y_pace = pd.read_csv(f'{indir}/MACLE2/label.csv')['censor_MACLE2']
test_model(X_pace, y_pace, outdir+'/evaluation/'+'MACLE2/')

########## PAD
X_pace = pd.read_csv(f'{indir}/PAD/normalized_counts.csv', index_col=0, header=0).T
y_pace = pd.read_csv(f'{indir}/PAD/label.csv')['PAD']
test_model(X_pace, y_pace, outdir+'/evaluation/'+'PAD/')