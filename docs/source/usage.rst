.. |PyPI| |Docs|
..
.. .. |PyPI| image:: https://img.shields.io/pypi/v/scanpy.svg
..    :target: https://pypi.org/project/scRFE/
.. .. |Docs| image:: https://readthedocs.com/projects/icb-scanpy/badge/?version=latest
..    :target: https://scRFE.readthedocs.io/en/latest/usage.html


Usage of scRFE
====================================
The key idea of scRFE is to find the features most important in describing the input data using one versus all classification.

1) **Read in AnnData Object and Filter/Normalize**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users will need to import Scanpy and read_h5ad from Scanpy for the following steps.

scRFE does not alter the AnnData object.
Thus, cells and genes must be filtered and normalized before running the code.
You may choose to filter out genes that appear in less than x number of cells, and cells that contain less than y number of genes.
Example code is below: (this will make the code run faster)


.. code:: bash

  sc.pp.filter_cells(AnnData, min_genes=250)
  sc.pp.filter_genes(AnnData, min_cells=3)


..

You may also want to normalize your data before running scRFE. Example code is below:

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
    scRFE (adata = AnnDataLiver, classOfInterest = 'age')


..


Default parameters are listed below, but users may manipulate them as they see fit.

    adata : anndata object
        The data file of interest.
    classOfInterest : str
        The class you will split the data by in the set of dataMatrix.obs.
    nEstimators : int, default = 5,000
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
