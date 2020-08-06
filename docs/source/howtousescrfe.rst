How To Use scRFE
====================================
The key idea of scRFE is to find the features most important in describing the input data using one versus all classification
by ranking feature importance.
Follow the steps below to run our transcription factor analysis listed in our scRFE preprint on bioRxiv.
To run analysis without transcription factors, simply omit Step 2.

Find a more detailed tutorial of how to use scRFE at https://github.com/czbiohub/scRFE/blob/master/scripts/scRFEtutorial.ipynb

1) **Read in AnnData Object and Filter/Normalize**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users will need to import Scanpy and read_h5ad from Scanpy for the following steps.

scRFE does not alter the AnnData object.
Thus, cells and genes must be filtered and normalized before running the code.
Users may choose to filter out genes that appear in less than x number of cells, and cells that contain less than y number of genes.
Example code is below: (this will make the code run faster)



.. code:: bash

  #imports
  from anndata import read_h5ad
  import pandas as pd
  from scRFE.scRFE import scRFE

  fileName = 'tabula-muris-senis-droplet-processed-official-annotations.h5ad' #AnnData object
  adata = read_h5ad(fileName)  # read in adata

  # basic filtering (optional)
  sc.pp.filter_cells(adata, min_genes=250)
  sc.pp.filter_genes(adata, min_cells=3)

..

Users may also want to normalize the data before running scRFE.
We did not notice significant differences between using raw
or normalized data, but users can evaulate this on their
own dataset. Example code is below:

.. code:: bash

  sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e5)
  sc.pp.log1p(adata)


..

Note that scRFE will automatically remove cells in the dataset if they don’t have a value for the class of interest.

2) **Read in Transcription Factor List and Subset With Adata**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash


  tfList = ‘GO_term_summary_20171110_222852.csv’ # transcription factor list
  mouseTFs = pd.read_csv(tfList) #read in TF list

  #subset adata to only have features from the transcription factor list
  adataTF = adata[:,adata.var_names[adata.var_names.isin(mouseTFs['Symbol'])]]


..



3) **Call scRFE**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Call scRFE specifying the AnnData object (adata) and the metadata column of adata
to split observations by (classOfInterest)

.. code:: bash

    # call scRFE to split observations by cell type
    scRFE(adata =  adataTF, classOfInterest = ‘cell_ontology_class’)


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
Users should use python’s pickle module to save the figure created.
Example code for running scRFE on the 3 month label only and plotting permutation importance is below.

.. code:: bash

  from scRFE.scRFE import makeOneForest
  from scRFE.scRFE import scRFEimplot
  3mForest = makeOneForest(dataMatrix=adata, classOfInterest=’age’, labelOfInterest=’3m’, nEstimators=1000,  randomState=0,  min_cells=15, keep_small_categories=True,   nJobs=-1, oobScore=True, Step=0.2, Cv=5, verbosity=True)
  fig = scRFEimplot(X_new = 3mForest[3], y = 3mForest[4])



..
