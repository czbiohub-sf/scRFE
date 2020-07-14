#!/usr/bin/env python
# coding: utf-8

# # scRFEimplot

# In[92]:


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


# In[65]:


def scRFEimplot(X_new,y):
    """
    Plots permutation importance of each feature selected by scRFE.
    Parameters
    ----------
    X_new : sparse matrix 
    Transformed array.
    y : pandas series
    Target labels.
    Returns
    -------
    plt : module matplotlib.pyplot
    Can be pickled, then saved as an image.
    """
    rf = RandomForestClassifier(random_state=0).fit(X_new, y)
    result = permutation_importance(rf, X_new.todense(), y, n_repeats=10, random_state=0,
        n_jobs=-1)
    fig, ax = plt.subplots()
    sorted_idx = result.importances_mean.argsort()
    ax.boxplot(result.importances[sorted_idx].T*100,
        vert=False, labels=range(X_new.shape[1]))
    ax.set_title("Permutation Importance of each feature")
    ax.set_ylabel("Features")
    fig.tight_layout()
    plt.savefig('plot.png', dpi=300, bbox_inches='tight') #trying to show

    plt.show()
    return plt


# In[66]:


# test3 = scRFEimplot(X_new = test1[3], y = test1[4])


# In[48]:


# type(test3)


# In[ ]:




