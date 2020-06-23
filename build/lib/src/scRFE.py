# Imports
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import read_h5ad
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import RFE
from sklearn.feature_selection import RFECV

adata = read_h5ad('/Users/madelinepark/Downloads/Limb_Muscle_facs.h5ad')
tiss = adata

def scRFE (tiss, feature='age', n_estimators=1000, random_state=0, n_jobs=-1, oob_score=True, test_size = 0.05, step=0.2, cv=5) :
    tiss.obs['feature_type_of_interest'] = 'rest'
    results_feature_cv = pd.DataFrame()
    for c in list(set(tiss.obs[feature])):
        print(c)
        clf = RandomForestClassifier(n_estimators=10, random_state=0, n_jobs=-1, oob_score=True)
        selector = RFECV(clf, step=0.2, cv=3, n_jobs=4) # step = % rounded down at each iteration
        feature_of_interest = c

        tiss.obs.loc[tiss.obs[tiss.obs[feature] == feature_of_interest].index,'feature_type_of_interest'] = feature_of_interest

        feat_labels = tiss.var_names
        X = tiss.X
        y = tiss.obs['feature_type_of_interest']

        print('training...')
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=0)

        clf.fit(X_train, y_train)
        selector.fit(X_train, y_train)
        feature_selected = feat_labels[selector.support_]

        print('result writing')

        column_headings = []
        column_headings.append(c)
        column_headings.append(c + '_gini')

        resaux = pd.DataFrame(columns=column_headings)
        resaux[c] = feature_selected
        resaux[c + '_gini'] = (selector.estimator_.feature_importances_)

        print(feature_selected)
        print (selector.estimator_.feature_importances_)

        results_feature_cv = pd.concat([results_feature_cv,resaux],axis=1)

        tiss.obs['feature_type_of_interest'] = 'rest'

    results_feature_cv

    scRFE(tiss=tiss, feature='age', n_estimators=10, step = 0.2, cv=2)
