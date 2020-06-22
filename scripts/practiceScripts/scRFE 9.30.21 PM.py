#!/usr/bin/env python
# coding: utf-8

# # scRFE

# In[ ]:





# In[84]:


# MENTION ONE VS ALL CLASSIFICATION in description


# In[117]:


# Imports 
import numpy as np
import pandas as pd
import scanpy as sc
import random
from anndata import read_h5ad
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import RFE
from sklearn.feature_selection import RFECV
import seaborn as sns
import matplotlib.pyplot as plt


# In[118]:


# read in anndata file 
adata = read_h5ad('/Users/madelinepark/Downloads/Kidney_facs.h5ad')


# In[119]:


def filterNormalize (dataMatrix, classOfInterest):
    np.random.seed(644685)
    sc.logging.print_versions()
    sc.settings.verbosity = 3      
    sc.logging.print_versions()
    tiss = dataMatrix
    tiss.obs['n_counts'] = tiss.X.sum(axis=1).A1
    sc.pp.filter_cells(tiss, min_genes=250)
    sc.pp.filter_genes(tiss, min_cells=3)
    tiss = tiss[tiss.obs['n_counts'] > 1500, :]
    sc.pp.normalize_per_cell(tiss, counts_per_cell_after=1e5)
    sc.pp.log1p(tiss)
    tiss.raw = tiss
    tiss = tiss[tiss.obs[classOfInterest]!='nan']
    return tiss


# In[120]:


# goal: get labels on a per class basis that will go into randomForest function for y
def getLabels (dataMatrix, classOfInterest): 
    """
    Gets labels on a per class basis that will inputted to the randomForest function
    
    Parameters
    ----------
    dataMatrix : anndata object
        The data file of interest
    classOfInterest : str
        The class you will split the data by in the set of dataMatrix.obs
    
    Returns
    -------
    labelsDict : dict
        Dictionary with labels for each class 
    """
    dataMatrix = filterNormalize (dataMatrix, classOfInterest)
    labelsDict = {}
    for label in np.unique(dataMatrix.obs[classOfInterest]):
        lists = []        
        for obs in dataMatrix.obs[classOfInterest]:
            if obs == label: 
                lists.append('A')
            else:
                lists.append('B')
        labelsDict[label] = lists #this is usually in line w if and else    
    return labelsDict


# In[121]:


def makeOneForest (dataMatrix, classOfInterest, labelOfInterest, nEstimators = 5000, 
                   randomState = 0,  nJobs = -1, oobScore = True, Step = 0.2, Cv = 5): 
    """
    Builds and runs a random forest for one label in a class of interest
    
    Parameters
    ----------
    dataMatrix : anndata object
        The data file of interest
    classOfInterest : str
        The class you will split the data by in the set of dataMatrix.obs
    labelOfInterest : str
        The specific label within the class that the random forezt will run a 
        "one vs all" classification on
    nEstimators : int
        The number of trees in the forest
    randomState : int
        Controls random number being used
    nJobs : int
        The number of jobs to run in parallel
    oobScore : bool
        Whether to use out-of-bag samples to estimate the generalization accuracy
    Step : float
        Corresponds to percentage of features to remove at each iteration
    Cv : int
        Determines the cross-validation splitting strategy
        
    Returns
    -------
    feature_selected : list
        list of top features from random forest
    selector.estimator_.feature_importances_ : list
        list of top ginis corresponding to to features
    
    """
    dataMatrix = filterNormalize (dataMatrix, classOfInterest)

    print('makeOneForest' + labelOfInterest)
    labelsDict = getLabels(dataMatrix, classOfInterest) 

    feat_labels = dataMatrix.var_names #this is equivalent of the genes
    X = dataMatrix.X
    y = labelsDict[labelOfInterest]
    print('Y')
    print(len(y))
    clf = RandomForestClassifier(n_estimators = nEstimators, random_state = randomState, 
                                 n_jobs = nJobs, oob_score = oobScore)
    selector = RFECV(clf, step = Step, cv = Cv)
    
    print('training...')
    clf.fit(X, y)
    selector.fit(X, y)
    feature_selected = feat_labels[selector.support_] 

    return feature_selected, selector.estimator_.feature_importances_ 


# In[122]:


def resultWrite (classOfInterest, results_df, labelOfInterest,
                feature_selected, feature_importance):
    print ('result writing')
    print(results_df)
    
    column_headings = [] 
    column_headings.append(labelOfInterest)
    column_headings.append(labelOfInterest + '_gini')
    resaux = pd.DataFrame(columns = column_headings)
    resaux[labelOfInterest] = feature_selected
    resaux[labelOfInterest + '_gini'] = feature_importance
    resaux = resaux.sort_values(by = [labelOfInterest + '_gini'], ascending = False)
    resaux.reset_index(drop = True, inplace = True)

    results_df = pd.concat([results_df, resaux], axis=1)
    return results_df 


# In[123]:


def scRFE(dataMatrix, classOfInterest, nEstimators = 5000, randomState = 0,  
                  nJobs = -1, oobScore = True, Step = 0.2, Cv = 5):
    """
    Builds and runs a random forest with one vs all classification for each label 
    for one class of interest
    
    Parameters
    ----------
    dataMatrix : anndata object
        The data file of interest
    classOfInterest : str
        The class you will split the data by in the set of dataMatrix.obs
    labelOfInterest : str
        The specific label within the class that the random forezt will run a 
        "one vs all" classification on
    nEstimators : int
        The number of trees in the forest
    randomState : int
        Controls random number being used
    nJobs : int
        The number of jobs to run in parallel
    oobScore : bool
        Whether to use out-of-bag samples to estimate the generalization accuracy
    Step : float
        Corresponds to percentage of features to remove at each iteration
    Cv : int
        Determines the cross-validation splitting strategy
        
    Returns
    -------
    results_df : pd.DataFrame
        Dataframe with results for each label in the class, formatted as 
        "label" for one column, then "label + gini" for the corresponding column
    
    """
    
    dataMatrix = filterNormalize (dataMatrix, classOfInterest)
    results_df = pd.DataFrame()
    for labelOfInterest in np.unique(dataMatrix.obs[classOfInterest]): #for timeliness    
        print( 'scRFE' + labelOfInterest)
        
        feature_selected, feature_importance = makeOneForest(dataMatrix, 
                                                             classOfInterest, 
                          labelOfInterest = labelOfInterest)
    
        results_df = resultWrite (classOfInterest, results_df, 
                            labelOfInterest = labelOfInterest, 
                    feature_selected = feature_selected,  
                    feature_importance = feature_importance)
        print(results_df.shape)
    return results_df


# # Run scRFE split by celltype

# In[125]:


kidney5000scRFECelltypeReIndex = scRFE (dataMatrix = adata, 
                                        classOfInterest = 'cell_ontology_class')

