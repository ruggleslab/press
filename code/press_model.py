###########################################################################
#
#                            press_model
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: press_model

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
experiment = "press_model_rfeada"
outdir = os.path.join(root_dir, "output", experiment) + "/"

# Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

## Set up some general items
utils.set_random_seed(420)
utils.hide_warnings()

# load the params json
params = pd.read_json('config/params.json', typ='series')

#======================== CODE ========================#
# Load the data
path = 'data/clean/'

genes = pd.read_table(params['geneset_path'][0]).values.flatten()

genes_df = pd.DataFrame(genes, columns=['Gene'])
genes_df.to_csv(outdir + 'genes.csv', index=False)

X_pace = pd.read_csv(path+'pace/features.csv', index_col=0, header=0)
y_pace = pd.read_csv(path+'pace/labels.csv').to_numpy()[:,0]

X_duke = pd.read_csv(path+'duke/features_group1.csv')
y_duke = pd.read_csv(path+'duke/labels_group1.csv').to_numpy()[:,0]

# X_test2 = pd.read_csv(path+'duke/features_group2.csv')
# y_test2 = pd.read_csv(path+'duke/labels_group2.csv').to_numpy()[:,0]

# subset to just the genes
X_pace = X_pace[genes]
X_duke = X_duke[genes]
# X_test2 = X_test2[genes]

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
    'rf0__max_features': ['log2', 'sqrt', 'auto'],
    'rf0__n_estimators': [10, 100, 200],
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
model = AdaBoostClassifier(learning_rate=0.1, n_estimators=1000)

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
#======================== Evaluation ========================#
######## Duke group 1
os.makedirs(outdir+'/evaluation/'+'duke', exist_ok=True)
duke_classification_report = metrics.classification_report(y_duke, model.predict(X_duke))
with open(outdir+'/evaluation/'+'duke/'+'classification_report.txt', 'w') as f:
    f.write(duke_classification_report)
plotting.plot_roc_curve_ci(model, X_duke, y_duke, title='Duke Cohort ROC Curve', save_path=outdir+'/evaluation/'+'duke/'+'duke_roc_curve.png')

########## HARP
# os.mkdir(outdir+'/evaluation/'+'harp')
