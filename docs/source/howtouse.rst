|PyPI| |Docs|

.. |PyPI| image:: https://img.shields.io/pypi/v/scanpy.svg
   :target: https://pypi.org/project/scRFE/
.. |Docs| image:: https://readthedocs.com/projects/icb-scanpy/badge/?version=latest
   :target: https://scRFE.readthedocs.io/en/latest/howtouse.html

How to use scRFE
=====================
The goal of scRFE is to find the genes that are most important in defining the transcriptome of the mouse for a given feature (age, celltype, etc). There are three main steps:

1) **Read in anndata file **
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use scanpy’s read_h5ad function to read in your anndata file.

This should be an anndata object contained in a data matrix X, annotations of observations (obs), variables (var), and unstructured annotations (uns).

2) **Set parameters**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N_estimators corresponds to the number of trees in your forest. We set the default to be 1000. More trees will take much longer to train, but give more accurate representations of the data. This random forest comes from sklearn, so please see sklearn.ensemble.RandomForestClassifier for more parameter details.

Step corresponds to the percent of your features that get taken out of your data set each time. The genes with the lowest ginis get taken out. This parameter comes from sklearn, so please see sklearn.feature_selection.RFECV for more parameter details.

Cv corresponds to the cross-validation splitting strategy which we set as a default to 5. The greater your cv, the more accurate your model becomes, however it will take much longer to turn. This parameter comes from sklearn, so please see sklearn.feature_selection.RFECV for more parameter details.

3) **Create for loop to split features**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We create the loop using the .loc function to ensure that we train the model separately for each value of our feature of interest. For example, we compare all 3 month data with itself, and the same for 24 month data, but not with each other when training the model.

After splitting the data so that the loop will loop through each value separately, we train the random forest classifier (clf) and our recursive feature selector (selector).

Finally, the loop finishes by writing the results to a data frame with one column for the gene name, and another for that gene’s corresponding gini.
