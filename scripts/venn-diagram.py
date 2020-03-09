#!/usr/bin/env python
# coding: utf-8

# # Visualization: Venn Diagram

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
from anndata import read_h5ad
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import RFE
from matplotlib import pyplot as plt


# In[2]:


import matplotlib_venn
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles


# # Venn Diagram

# In[3]:


# Read in SORTED results for each feature (ex: age)
# Highest gini scores and their corresponding genes at top, fill the column in descending order
sorted_24 = pd.read_csv('/Users/madelinepark/src2/maca-data-analysis/rf-rfe-results/cv_24m_facs_sorted.csv')
sorted_3 = pd.read_csv('/Users/madelinepark/src2/maca-data-analysis/rf-rfe-results/cv_3m_facs_sorted.csv')


# In[4]:


# compare 3m vs 24m results
compare_3_24 = venn2([set(sorted_24['24m']), set(sorted_3['3m'])], set_labels=('3m','24m'))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




