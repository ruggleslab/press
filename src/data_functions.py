import pandas as pd
import sys
import os

# Function to load the data
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

# Function to normalize the data
def normalize_data(X, scaler):
    """
    Normalize the input data.

    Parameters
    ----------
    X : pd.DataFrame
        The input data to normalize.
    scaler : object
        The scaler to use for normalization.

    Returns
    -------
    X : pd.DataFrame
        The normalized input data.
    """
    X = pd.DataFrame(scaler.transform(X), columns=X.columns, index=X.index)
    return X