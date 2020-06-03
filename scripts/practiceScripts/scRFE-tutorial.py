#!/usr/bin/env python
# coding: utf-8

# # scRFE Tutorial
# 

# Here we present an example of how to use scRFE. We analyze the Limb Muscle Facs data from the Tabula-Muris-Senis dataset that is available on Figshare. We split the data by age.

# More features were selected than ideal in this model, because we used a very small number of estimators and a low CV score, for time's sake. This results are not accurate though, and we recommend running the code with 1000 estimators and CV>=5 with an EC2 instance.

# ### Imports 

# In[2]:


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


# ### Read in anndata file 

# In[3]:


adata = read_h5ad('/Users/madelinepark/Downloads/Limb_Muscle_facs.h5ad')
tiss = adata


# In[4]:


list(set(tiss.obs['age']))


# ### Run scRFE

# we decreased n_estimators and cv so that the code will run faster, but you should increase both before using

# In[5]:


tiss.obs['age_type_of_interest'] = 'rest'
results_age_cv = pd.DataFrame()

for c in list(set(tiss.obs['age'])): 
    print(c)
    clf = RandomForestClassifier(n_estimators=10, random_state=0, n_jobs=-1, oob_score=True)
    selector = RFECV(clf, step=0.2, cv=3, n_jobs=4) # step = % rounded down at each iteration  
    age_of_interest = c
    
    tiss.obs.loc[tiss.obs[tiss.obs['age'] == age_of_interest].index,'age_type_of_interest'] = age_of_interest
    
    feat_labels = tiss.var_names 
    X = tiss.X
    y = tiss.obs['age_type_of_interest']
    
    print('training...')
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=0) 
    
    clf.fit(X_train, y_train)
    selector.fit(X_train, y_train)
    print(type(feat_labels)) #adding this to test
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
    
    results_age_cv = pd.concat([results_age_cv,resaux],axis=1)
    
    tiss.obs['age_type_of_interest'] = 'rest'
    
results_age_cv


# In[ ]:





# In[ ]:





# In[ ]:




