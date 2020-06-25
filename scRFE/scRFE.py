#!/usr/bin/env python
# coding: utf-8

# # scRFE

# In[3]:


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
import scanpy.external as sce
import logging as logg


# In[84]:


def columnToString (dataMatrix):
    cat_columns = dataMatrix.obs.select_dtypes(['category']).columns
    dataMatrix.obs[cat_columns] = dataMatrix.obs[cat_columns].astype(str)

    return dataMatrix


# In[101]:


def filterNormalize (dataMatrix, classOfInterest):
    np.random.seed(644685)
#     sc.pp.filter_cells(dataMatrix, min_genes=0)
#     sc.pp.filter_genes(dataMatrix, min_cells=0)
    dataMatrix = dataMatrix[dataMatrix.obs[classOfInterest]!='nan']
    dataMatrix = dataMatrix[~dataMatrix.obs[classOfInterest].isna()]
    print ('na data removed')
    return dataMatrix


# In[86]:


def labelSplit (dataMatrix, classOfInterest, labelOfInterest):
    dataMatrix = filterNormalize (dataMatrix, classOfInterest)
    dataMatrix.obs['classification_group'] = 'B'
    dataMatrix.obs.loc[dataMatrix.obs[dataMatrix.obs[classOfInterest]==labelOfInterest]
                   .index,'classification_group'] = 'A' #make labels based on A/B of
#     classofInterest
    return dataMatrix


# In[92]:


def downsampleToSmallestCategory(dataMatrix, random_state, min_cells,
                                 keep_small_categories,
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


# In[93]:


def makeOneForest (dataMatrix, classOfInterest, labelOfInterest, nEstimators,
                   randomState,  min_cells, keep_small_categories,
                   nJobs, oobScore, Step, Cv):

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
    splitDataMatrix = labelSplit (dataMatrix, classOfInterest, labelOfInterest)

    downsampledMatrix = downsampleToSmallestCategory (dataMatrix = splitDataMatrix,
    random_state = randomState, min_cells = min_cells,
        keep_small_categories = keep_small_categories,
        classOfInterest = 'classification_group', )

    feat_labels = downsampledMatrix.var_names
    X = downsampledMatrix.X
    y = downsampledMatrix.obs['classification_group'] #'A' or 'B' labels from labelSplit

    clf = RandomForestClassifier(n_estimators = nEstimators, random_state = randomState,
                                 n_jobs = nJobs, oob_score = oobScore)

    selector = RFECV(clf, step = Step, cv = Cv)

    clf.fit(X, y)
    selector.fit(X, y)
    feature_selected = feat_labels[selector.support_]
    dataMatrix.obs['classification_group'] = 'B'

    return feature_selected, selector.estimator_.feature_importances_


# In[94]:


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


# In[99]:


def scRFE (adata, classOfInterest, nEstimators = 5000, randomState = 0, min_cells = 15,
        keep_small_categories = True, nJobs = -1, oobScore = True, Step = 0.2, Cv = 5):

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

    Returns
    -------
    results_df : pd.DataFrame
        Dataframe with results for each label in the class, formatted as
        "label" for one column, then "label + gini" for the corresponding column

    """
    dataMatrix = adata.copy()
    dataMatrix = columnToString (dataMatrix)
    dataMatrix = filterNormalize (dataMatrix, classOfInterest)
    results_df = pd.DataFrame()

    for labelOfInterest in sorted(np.unique(dataMatrix.obs[classOfInterest])):
        dataMatrix_labelOfInterest = dataMatrix.copy()

        feature_selected, feature_importance =  makeOneForest(
            dataMatrix = dataMatrix_labelOfInterest, classOfInterest = classOfInterest,
            labelOfInterest = labelOfInterest,
            nEstimators = nEstimators, randomState = randomState,  min_cells = min_cells,
            keep_small_categories = keep_small_categories, nJobs = nJobs,
            oobScore = oobScore, Step = Step, Cv = Cv)

        results_df = resultWrite (classOfInterest, results_df,
                            labelOfInterest = labelOfInterest,
                    feature_selected = feature_selected,
                    feature_importance = feature_importance)


    return results_df
