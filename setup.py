from setuptools import setup, find_packages

setup(
    name = "scRFE",
    version = "0.1",
    keywords = ("pip", "single cell", "scRFE"),
    description = "Single-cell identity definition using random forest modelling and recursive feature elimination",
    license = "MIT Licence",

    url = "https://github.com/wangshenguiuc/OnClass",
    author = "Madeline Park",
    author_email = "madeline.park@czbiohub.org",

    packages = find_packages(),
    include_package_data = True,
    platforms = "any",
    install_requires=[
		'anndata>=0.6.22.post1',
        'umap-learn>=0.3.10',
        'matplotlib>=2.0.2',
        'numpy>=1.16.4',
        'scikit-learn>=0.21.3'

    ]
)
