How To Use scRFE
====================================
The key idea of scRFE is to find the features most important in describing the input data using one versus all classification.

1) **Read in AnnData Object and Filter/Normalize**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users will need to import Scanpy and read_h5ad from Scanpy for the following steps.

scRFE does not alter the AnnData object.
Thus, cells and genes must be filtered and normalized before running the code.
Users may choose to filter out genes that appear in less than x number of cells, and cells that contain less than y number of genes.
Example code is below: (this will make the code run faster)



.. code:: bash

  sc.pp.filter_cells(AnnData, min_genes=250)
  sc.pp.filter_genes(AnnData, min_cells=3)


..

Users may also want to normalize your data before running scRFE. Example code is below:

.. code:: bash

  sc.pp.normalize_per_cell(AnnData, counts_per_cell_after=1e5)
  sc.pp.log1p(AnnData)


..

Note that scRFE will automatically remove cells in the dataset if they don’t have a value for the class of interest.

2) **Call scRFE and Return Results DataFrame**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once the data is ready, call scRFE, specifying the filtered and normalized
AnnData object and the classOfInterest that the data will be classified by.
Example code is below:

.. code:: bash

    from scRFE.scRFE import scRFE
    scRFEdf = scRFE (adata = AnnDataLiver, classOfInterest = 'age')


..


Default parameters are listed below, but users may manipulate them as they see fit.

    adata : anndata object
        The input .h5ad file of interest.
    classOfInterest : str
        The class you will classify the data by in the set of adata.obs.
    nEstimators : int, default = 1,000
        The number of trees in the forest.
    randomState : int, default = 0
        Controls random number being used.
    min_cells : int, default = 15
        Minimum number of cells in a given class to downsample.
    keep_small_categories : bool, default = True
        Whether to keep classes with small number of observations, or to remove.
    nJobs : int, default = -1
        The number of jobs to run in parallel.
    oobScore : bool, default = True
        Whether to use out-of-bag samples to estimate the generalization accuracy.
    Step : float, default = 0.2
        Corresponds to percentage of features to remove at each iteration.
    Cv : int, default = 5
        Determines the k-fold cross-validation splitting strategy.
    verbosity : bool, default = True
        Whether to include print statements.

3) **Plot permutation importance for each feature.**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If users are interested in seeing each feature’s importance for a given label within the classOfInterest, they can run makeOneForest to extract X_new and y.
X_new is a sparse matrix containing the selected features for one label in the classOfInterest.
y is a pandas series with the target labels. These are the two inputs to scRFEimplot, which will then plot the importance of each feature.
Users should use python’s pickle module to save the figure created. Example code is below:

.. code:: bash

  from scRFE.scRFE import makeOneForest
  from scRFE.scRFE import scRFEimplot
  3mForest = makeOneForest(dataMatrix=adataLiver, classOfInterest=’age’, labelOfInterest=’3m’, nEstimators=1000,  randomState=0,  min_cells=15, keep_small_categories=True,   nJobs=-1, oobScore=True, Step=0.2, Cv=5, verbosity=True)
  fig = scRFEimplot(X_new = 3mForest[3], y = 3mForest[4])



..
