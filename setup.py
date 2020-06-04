from setuptools import setup, find_packages

setup(
    name = "scRFE",
    version = "0.6",
    keywords = ("pip", "single cell", "scRFE"),
    description = "Single-cell identity definition using random forest modelling and recursive feature elimination",
    long_description = "Single-cell identity one vs all classification using random forest modelling and recursive feature elimination",
    license = "MIT Licence",

    url = "https://github.com/czbiohub/scRFE",
    author = "Madeline Park",
    author_email = "madeline.park@czbiohub.org",

    packages = find_packages(),
    include_package_data = True,
    platforms = "any",
    install_requires=[
		'anndata>= 0.6.21',
        'matplotlib>=3.1.1',
        'numpy>=1.16.4',
        'scikit-learn>=0.20.3',
        'scanpy>=1.4.3',
        'pandas>=0.24.2'

    ]
)
