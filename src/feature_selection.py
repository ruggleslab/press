import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.feature_selection import RFECV, SequentialFeatureSelector
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

class RFESelector:
    def __init__(self, model, scoring="roc_auc"):
        """
        Initialize the RFESelector class.

        Parameters:
        model (object): The model to use for feature selection.
        scoring (str): Scoring metric for feature selection. Default is "roc_auc".
        """
        self.model = model
        self.scoring = scoring
        self.rfecv = []
        self.X = None
        self.y = None

    def fit(self, X, y):
        """
        Perform feature selection using Recursive Feature Elimination (RFE).

        Parameters:
        X (DataFrame): The feature matrix.
        y (Series): The target variable.

        Returns:
        None
        """
        # Create recursive feature elimination object
        cv = StratifiedKFold(5)  # Cross-validation generator
        self.rfecv = RFECV(estimator=self.model, step=1, cv=cv, scoring=self.scoring)

        # Fit recursive feature elimination object
        self.rfecv.fit(X, y)
        self.X = X
        self.y = y

        print(f"Optimal number of features: {self.rfecv.n_features_}")
        
    def plot_features(self, output_dir=None):
        """
        Plot the number of features vs. cross-validation scores.

        Returns:
        None
        """
        # Check if get_features has been called
        if self.rfecv is None:
            raise RuntimeError("You must call get_features before calling plot_features.")

        os.makedirs(output_dir, exist_ok=True)

        # Plot number of features VS. cross-validation scores
        n_scores = len(self.rfecv.cv_results_["mean_test_score"])
        plt.figure()
        plt.plot(
            range(1, n_scores + 1),
            self.rfecv.cv_results_["mean_test_score"],
            color="b", label="mean", lw=1,
        )
        plt.fill_between(
            range(1, n_scores + 1),
            self.rfecv.cv_results_["mean_test_score"] - self.rfecv.cv_results_["std_test_score"],
            self.rfecv.cv_results_["mean_test_score"] + self.rfecv.cv_results_["std_test_score"],
            color="r", alpha=0.2, label="std",
        )
        plt.legend()
        if output_dir:
            plt.savefig(os.path.join(output_dir, "rfe_features.png"))
        else:
            plt.show()

    def save_results(self, output_dir):
        """
        Save the selected features to an output directory.

        Parameters:
        output_dir (str): Path to the output directory.

        Returns:
        None
        """
        # Check if get_features has been called
        if self.rfecv is None:
            raise RuntimeError("You must call get_features before calling save_results.")

        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Get the selected features
        selected_features = self.X.columns[self.rfecv.support_]

        # Save the selected features to a CSV file
        selected_features.to_series().to_csv(os.path.join(output_dir, 'selected_features.csv'), index=False)

        print(f'Saved {len(selected_features)} selected features to {output_dir}/selected_features.csv.')

class SFSSelector:
    def __init__(self, model, k_features=3, scoring='accuracy'):
        """
        Initialize the SFSSelector class.

        Parameters:
        model (object): The model to use for feature selection.
        k_features (int): The number of features to select. Default is 3.
        scoring (str): Scoring metric for feature selection. Default is 'accuracy'.
        """
        self.model = model
        self.k_features = k_features
        self.scoring = scoring
        self.sfs = None
        self.X = None
        self.y = None

    def fit(self, X, y):
        """
        Perform feature selection using Sequential Feature Selection (SFS).

        Parameters:
        X (DataFrame): The feature matrix.
        y (Series): The target variable.

        Returns:
        None
        """
        self.sfs = SequentialFeatureSelector(self.model, n_features_to_select=self.k_features, scoring=self.scoring)
        self.X = X
        self.y = y
        self.sfs.fit(X, y)

        print(f"Selected features: {X.columns[self.sfs.get_support()].tolist()}")

    def plot_results(self, output_dir=None):
        """
        Plot the results of the Sequential Feature Selection.

        Returns:
        None
        """
        # Check if fit has been called
        if self.sfs is None:
            raise RuntimeError("You must call fit before calling plot_results.")
        
        os.makedirs(output_dir, exist_ok=True)
        plt.plot(
            range(1, len(self.sfs.get_support(indices=True)) + 1), 
            cross_val_score(self.model, self.X[:, self.sfs.get_support(indices=True)], self.y, scoring=self.scoring, cv=StratifiedKFold(5))
            )
        plt.title('Sequential Forward Selection')
        plt.xlabel('Number of features selected')
        plt.ylabel('Cross validation score')
        plt.grid()
        if output_dir:
            plt.savefig(os.path.join(output_dir, "sfs_features.png"))
        else:
            plt.show()

    def save_results(self, output_dir):
        """
        Save the selected features to an output directory.

        Parameters:
        output_dir (str): Path to the output directory.

        Returns:
        None
        """
        # Check if fit has been called
        if self.sfs is None:
            raise RuntimeError("You must call fit before calling save_results.")

        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Save the selected features to a CSV file
        pd.Series(self.X.columns[self.sfs.get_support()].tolist()).to_csv(os.path.join(output_dir, 'selected_features.csv'), index=False)

        print(f'Saved {len(self.X.columns[self.sfs.get_support()].tolist())} selected features to {output_dir}/selected_features.csv.')
        
# Millipede Selector ### NOT CURRENTLY FINISHED ###
# class MillipedeSelector:
#     def __init__(self, T=2000, T_burnin=1000):
#         """
#         Initialize the MillipedeSelector class.
        
#         Parameters:
#         T (int): The number of iterations. Default is 2000.
#         T_burnin (int): The number of burn-in iterations. Default is 1000.
#         """
#         self.millipede = None
#         self.T = T
#         self.T_burnin = T_burnin
#         self.results = None
#         self.dat = None

#     def fit(self, X, y):
#         """
#         Perform feature selection using Millipede.

#         Parameters:
#         X (DataFrame): The feature matrix.
#         y (Series): The target variable.

#         Returns:
#         None
#         """
#         self.dat = pd.concat([X, y], axis=1)
#         self.dat = self.dat.rename(columns={0: 'label'})
#         self.millipede = BernoulliLikelihoodVariableSelector(self.dat, 'label')
#         self.millipede.run(T=self.T, T_burnin=self.T_burnin)
#         self.results = self.millipede.summary()

#         print(f"Selected features: {X.columns[self.millipede.get_support()].tolist()}")

#     def plot_results(self, output_dir=None):
#         """
#         Plot the results of the Millipede feature selection.

#         Returns:
#         None
#         """
#         # Check if fit has been called
#         if self.millipede is None:
#             raise RuntimeError("You must call fit before calling plot_results.")

#         plt.plot(
#             range(1, len(self.millipede.get_support(indices=True)) + 1), 
#             cross_val_score(self.model, self.X[:, self.millipede.get_support(indices=True)], self.y, scoring=self.scoring, cv=StratifiedKFold(5))
#             )
#         plt.title('Millipede Feature Selection')
#         plt.xlabel('Number of features selected')
#         plt.ylabel('Cross validation score')
#         plt.grid()
#         if output_dir:
#             plt.savefig(os.path.join(output_dir, "millipede_features.png"))
#         else:
#             plt.show()

#     def save_results(self, output_dir):
#         """
#         Save the selected features to an output directory.

#         Parameters:
#         output_dir (str): Path to the output directory.

#         Returns:
#         None
#         """
#         # Check if fit has been called
#         if self.millipede is None:
#             raise RuntimeError("You must call fit before calling save_results.")
#         summ = self.millipede.summary
#         # Create the output directory if it doesn't exist
#         os.makedirs(output_dir, exist_ok=True)
#         summ.to_csv(os.path.join(output_dir, 'selected_features.csv'), index=False)