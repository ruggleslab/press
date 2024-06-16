# Script to run the press model on a dataset
import argparse
import pandas as pd
from sklearn import preprocessing
from joblib import load
from MattTools import utils

def predict(data_path, model):
    data = pd.read_csv(data_path, index_col=0)
    print("Data shape: ", data.shape)
    print(data.head())
    # scaling the data
    scaler = preprocessing.StandardScaler()
    data = pd.DataFrame(
        scaler.fit_transform(data), 
        index=data.index, 
        columns=data.columns
        )
    
    # get the predictions and scores
    preds = model.predict_proba(data)[::,1]
    # if the shape of the data is [n, < 3] then don't scale
    if data.shape[1] < 3:
        print("Using predefined scaling.")
        # predefined mean is 0.178 and std is 0.858
        scores = (preds - 0.178) / 0.858
    else:
        scores = preprocessing.scale(preds)

    # make a dataframe with the scores and preds
    out = pd.DataFrame(
        {   
            "preds": preds,
            "scores": scores
        },
        index=data.index
    )
    return out

def check_press(data_path):
    data = pd.read_csv(data_path, index_col=0)
    # return an error if the data is not the right shape
    if data.shape[1] != 451:
        raise ValueError("Data must have 451 columns (genes).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to predict using the press model.")
    parser.add_argument("--data", required=True, help="Path to the data file in the format of Samples x Genes (### x 451)")
    parser.add_argument("--out", required=True, help="Path to the output file")
    args = parser.parse_args()
    
    utils.hide_warnings()
    model = load("models/jobs/press451.joblib")
    print("Checking data...")
    check_press(args.data)
    
    print("Predicting...")
    pred_df = predict(args.data, model)
    print("Done predicting.")

    print("Writing predictions to: ", args.out)
    pred_df.to_csv(args.out)
