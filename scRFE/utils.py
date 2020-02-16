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

def loc_split(tiss, feature='age'):
    tiss.obs['feature_type_of_interest'] = 'rest'
    results_feature_cv = pd.DataFrame()
    for c in list(set(tiss.obs[feature])):
        feature_of_interest = c
        tiss.obs.loc[tiss.obs[tiss.obs[feature] == feature_of_interest].index,'feature_type_of_interest'] = feature_of_interest
        feat_labels = tiss.var_names
        X = tiss.X
        y = tiss.obs['feature_type_of_interest']
    return X, y, feature, feat_labels

def train(X, y, test_size, random_state):
    print('training...')
    feat_labels = loc_split(tiss=tiss, feature='age')[3]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=0)
    clf = RandomForestClassifier(n_estimators=1000, random_state=0, n_jobs=-1, oob_score=True)
    selector = RFECV(clf, step=0.2, cv=3, n_jobs=4)
    clf.fit(X_train, y_train)
    selector.fit(X_train, y_train)
    feature_selected = feat_labels[selector.support_]
    return selector, clf, feat_labels, feature_selected, selector.estimator_.feature_importances_, X_train, X_test, y_train, y_test

def scRFE (tiss=tiss, X=tiss.X, feature='age', n_estimators=1000, random_state=0, n_jobs=-1, oob_score=True, test_size = 0.05, step=0.2, cv=5) :
    for c in list(set(tiss.obs[feature])):
        print(c)
        X = loc_split(tiss=tiss, feature= feature)[0] #change age to feature
        y = loc_split(tiss=tiss, feature= feature)[1]
        feature = loc_split(tiss=tiss, feature= feature)[2]
        feat_labels = loc_split(tiss=tiss, feature= feature)[3]
        feature_selected = train(X, y, test_size, random_state)[1]
        X_train = train(X, y, test_size, random_state)[2]
        X_test = train(X, y, test_size, random_state)[3]
        y_train = train(X, y, test_size, random_state)[4]
        y_test = train(X, y, test_size, random_state)[5]
        result_write(c, feature_selected, tiss, feature=feature,test_size=0.05, random_state=0)
        tiss.obs['age_type_of_interest'] = 'rest'
#         not sorted

scRFE(tiss=tiss, feature='age', n_estimators=1000, random_state=0, n_jobs=-1, oob_score=True, test_size = 0.05, step=0.2, cv=5)
