#!/usr/bin/env python
# coding: utf-8

# In[168]:


# import dependencies
import numpy as np
import pandas as pd
import scanpy as sc
import random
import logging as logg
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.metrics import accuracy_score
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt
from tqdm import tqdm 


# In[146]:


# transform all category columns in string columns
def columnToString (dataMatrix):
    cat_columns = dataMatrix.obs.select_dtypes(['category']).columns
    dataMatrix.obs[cat_columns] = dataMatrix.obs[cat_columns].astype(str)
    return dataMatrix


# In[147]:


# remove observations that are NaN for the category
def filterNormalize (dataMatrix, classOfInterest, verbosity):
    np.random.seed(644685)
    dataMatrix = dataMatrix[dataMatrix.obs[classOfInterest]!='nan']
    dataMatrix = dataMatrix[~dataMatrix.obs[classOfInterest].isna()]
    if verbosity == True:
        print ('Removed NaN observations in the selected category')
    return dataMatrix


# In[148]:


# set the A/B labels for classification
def labelSplit (dataMatrix, classOfInterest, labelOfInterest, verbosity):
    dataMatrix = filterNormalize (dataMatrix, classOfInterest, verbosity)
    dataMatrix.obs['classification_group'] = 'B'
    dataMatrix.obs.loc[dataMatrix.obs[dataMatrix.obs[classOfInterest]==labelOfInterest]
                   .index,'classification_group'] = 'A' #make labels based on A/B of classofInterest
    return dataMatrix


# In[149]:


# downsample observations to balance the groups
def downsampleToSmallestCategory(dataMatrix, random_state, min_cells,
                                 keep_small_categories, verbosity,
                                 classOfInterest = 'classification_group',
) -> sc.AnnData:
    """
    returns an annData object in which all categories in 'classOfInterest' have
    the same size
    classOfInterest
        column with the categories to downsample
    min_cells
        Minimum number of cells to downsample.
        Categories having less than `min_cells` are discarded unless
        keep_small_categories is True
    keep_small_categories
        Be default categories with less than min_cells are discarded.
        Set to true to keep them
    """
    counts = dataMatrix.obs[classOfInterest].value_counts(sort=False)
    if len(counts[counts < min_cells]) > 0 and keep_small_categories is False:
        logg.warning(
            "The following categories have less than {} cells and will be "
            "ignored: {}".format(min_cells, dict(counts[counts < min_cells]))
        )
    min_size = min(counts[counts >= min_cells])
    sample_selection = None
    for sample, num_cells in counts.items():
        if num_cells <= min_cells:
            if keep_small_categories:
                sel = dataMatrix.obs.index.isin(
                    dataMatrix.obs[dataMatrix.obs[classOfInterest] == sample].index)
            else:
                continue
        else:
            sel = dataMatrix.obs.index.isin(
                dataMatrix.obs[dataMatrix.obs[classOfInterest] == sample]
                .sample(min_size, random_state=random_state)
                .index
            )
        if sample_selection is None:
            sample_selection = sel
        else:
            sample_selection |= sel
    logg.info(
        "The cells in category {!r} had been down-sampled to have each {} cells. "
        "The original counts where {}".format(classOfInterest, min_size, dict(counts))
    )
    return dataMatrix[sample_selection].copy()


# In[155]:


# build the random forest classifier and perform feature elimination
def makeOneForest (dataMatrix, classOfInterest, labelOfInterest, nEstimators,
                   randomState,  min_cells, keep_small_categories,
                   nJobs, oobScore, Step, Cv, verbosity):
    #need to add verbose arg details 
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
    verbosity : bool 
        Whether to include print statements.
    Returns
    -------
    feature_selected : list
        list of top features from random forest
    selector.estimator_.feature_importances_ : list
        list of top ginis corresponding to to features
    """
    splitDataMatrix = labelSplit (dataMatrix, classOfInterest, labelOfInterest, verbosity)

    downsampledMatrix = downsampleToSmallestCategory (dataMatrix = splitDataMatrix, 
    random_state = randomState, min_cells = min_cells, 
        keep_small_categories = keep_small_categories, verbosity = verbosity,
        classOfInterest = 'classification_group' )
    
    if verbosity == True:
        print(labelOfInterest)
        print(pd.DataFrame(downsampledMatrix.obs.groupby(['classification_group',classOfInterest])[classOfInterest].count()))

    feat_labels = downsampledMatrix.var_names
    X = downsampledMatrix.X
    y = downsampledMatrix.obs['classification_group'] #'A' or 'B' labels from labelSplit

    clf = RandomForestClassifier(n_estimators = nEstimators, random_state = randomState,
                                 n_jobs = nJobs, oob_score = oobScore)

    Cv = StratifiedKFold(Cv)
    selector = RFECV(clf, step = Step, cv = Cv, scoring='f1_weighted', min_features_to_select=2)

    clf.fit(X, y)
    selector.fit(X, y)
    feature_selected = feat_labels[selector.support_]
    dataMatrix.obs['classification_group'] = 'B'
    
    X_new = selector.fit_transform(X, y)
    selector.fit(X_new, y) 
    score = selector.score(X_new, y)
    feature_selected = feature_selected[selector.support_]

    return feature_selected, selector.estimator_.feature_importances_,score,X_new,y


# In[156]:


# write the results
def resultWrite (classOfInterest, results_df, labelOfInterest,
                feature_selected, feature_importance):

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


# In[164]:


# main scRFE function
def scRFE (adata, classOfInterest, nEstimators = 5000, randomState = 0, min_cells = 15,
        keep_small_categories = True, nJobs = -1, oobScore = True, Step = 0.2, Cv = 5, 
          verbosity = True):    
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
    min_cells : int
        Minimum number of cells in a given class to downsample.
    keep_small_categories : bool
        Whether to keep classes with small number of observations, or to remove.
    nJobs : int
        The number of jobs to run in parallel
    oobScore : bool
        Whether to use out-of-bag samples to estimate the generalization accuracy
    Step : float
        Corresponds to percentage of features to remove at each iteration
    Cv : int
        Determines the cross-validation splitting strategy
    verbosity : bool 
        Whether to include print statements.
    Returns
    -------
    results_df : pd.DataFrame
        Dataframe with results for each label in the class, formatted as
        "label" for one column, then "label + gini" for the corresponding column.
    score_df: dict
        Score for each label in classOfInterest.
    """
    
    dataMatrix = adata.copy()
    dataMatrix = columnToString (dataMatrix)
    dataMatrix = filterNormalize (dataMatrix, classOfInterest, verbosity)
    results_df = pd.DataFrame()
    
    score_df = {}
    for i in tqdm(range(len(adata.obs[classOfInterest]))[0:3]):

        for labelOfInterest in np.unique(dataMatrix.obs[classOfInterest])[0:3]:
            
            dataMatrix_labelOfInterest = dataMatrix.copy()

            feature_selected, feature_importance, model_score, X_new, y =  makeOneForest(dataMatrix = dataMatrix, 
                classOfInterest = classOfInterest, labelOfInterest = labelOfInterest, 
                nEstimators = nEstimators, randomState = randomState,  min_cells = min_cells, 
                    keep_small_categories = keep_small_categories,
                       nJobs = nJobs, oobScore = oobScore, Step= Step, Cv=Cv, verbosity=verbosity)

            results_df = resultWrite (classOfInterest, results_df,
                                labelOfInterest = labelOfInterest,
                        feature_selected = feature_selected,
                        feature_importance = feature_importance)

            score_df[labelOfInterest] = model_score
        pass

    return results_df,score_df

