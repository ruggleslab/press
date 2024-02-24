# This script is mainly to prep validation data for use with a given gene signature. It will load the data, normalize it, and then save it to a file for use in the validation step.

# +
import os
import pandas as pd
from joblib import load, dump
from sklearn.preprocessing import StandardScaler
import sys

# Function to load the data from a given counts and metadata file
def load_data(counts, labels):
    """
    Load the data from the input files.

    Parameters
    ----------
    counts : str
        Path to the input counts.
    labels : str
        Path to the input labels.

    Returns
    -------
    X : pd.DataFrame
        The input counts.
    y : pd.DataFrame
        The input labels.
    """
    X = pd.read_csv(counts, index_col=0)
    y = pd.read_csv(labels)
    return X, y