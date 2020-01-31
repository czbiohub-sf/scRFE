|PyPI| |Docs|

.. |PyPI| image:: https://img.shields.io/pypi/v/scanpy.svg
   :target: https://pypi.org/project/OnClass/
.. |Docs| image:: https://readthedocs.com/projects/icb-scanpy/badge/?version=latest
   :target: https://onclass.readthedocs.io/en/latest/introduction.html

API
=====================
**Main scRFEE **
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:param input_set: anndata object, .h5ad file from figshare

 :param feature: str, pick any column from tiss.obs (ex: age, cell_ontology_class)

 :param n_estimators : integer, optional (default=1000), the number of trees in the forest.
     see sklearn.ensemble.RandomForestClassifier for parameter details.

 :random_state : int, (default=0), see sklearn.model_selection.train_test_split for parameter details.

 :param n_jobs : int or None, (default=-1), she number of jobs to run in parallel for both fit and predict. -1 means using all processors. see sklearn.ensemble.RandomForestClassifier for parameter details.

 :param oob_score : bool (default=True), uses out-of-bag samples to estimate the generalization accuracy. see sklearn.ensemble.RandomForestClassifier for parameter details.

 :param step: int or float, optional (default=0.2), step corresponds to the percentage of features to remove at each iteration. see sklearn.feature_selection.RFECV for parameter details.

 :param cv: int, cross-validation generator or an iterable, optional (default = 5), determines the cross-validation splitting strategy. greater cv more accurate, takes longer to run. see sklearn.feature_selection.RFECV for parameter details.

 :param test_size: float, int or None (default = 0.05), represents the proportion of the dataset to include in the test split. see sklearn.model_selection.train_test_split for parameter details.
