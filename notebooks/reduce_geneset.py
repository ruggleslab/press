import argparse
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel

# Create an argument parser
parser = argparse.ArgumentParser(description='Perform feature selection.')
parser.add_argument('--counts', type=str, help='path to the input counts')
parser.add_argument('--labels', type=str, help='path to the input labels')

# Parse the arguments
args = parser.parse_args()

# Load the data
X = pd.read_csv(args.dataset)
y = pd.read_csv(args.labels)

# Define the model
model = RandomForestClassifier(n_estimators=100)

# Perform feature selection
selector = SelectFromModel(estimator=model)
X_new = selector.fit_transform(X, y)

# Get the selected feature indices
selected_indices = selector.get_support(indices=True)
selected_features = X.columns[selected_indices]

# Save the selected features
selected_features.to_csv('selected_features.csv', index=False)
X_new.to_csv('selected_features.csv', index=True)