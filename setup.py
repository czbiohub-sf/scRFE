from setuptools import setup, find_packages
setup(
  name = 'scRFE',         # How you named your package folder (MyLib)
  packages = ['scRFE'],   # Chose the same as "name"
  version = '1.4.1',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = "Single-cell identity definition using one vs all random forest classification and recursive feature elimination",
  license = "MIT Licence",
  author = 'Madeline Park',                   # Type in your name
  author_email = 'madeline.park@czbiohub.org',      # Type in your E-Mail
  url = "https://github.com/czbiohub/scRFE",   # Provide either the link to your github or to your website
  download_url = 'https://github.com/czbiohub/scRFE/archive/1.4.1.tar.gz',    # I explain this later on
  keywords = ("pip", "single cell", "scRFE"),   # Keywords that define your package best

  packages = find_packages(),
  include_package_data = True,
  platforms = "any",
  install_requires=[
    'anndata>= 0.6.21',
    'matplotlib>=3.1.1',
    'numpy>=1.16.4',
    'scikit-learn>=0.20.3',
    'scanpy>=1.4.3',
    'pandas>=0.24.2',
    'seaborn>=0.9.0',
    ],

  classifiers=[
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3.6',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.7',
  ],
)
#
# setup(
#     name = "scRFE",
#     version = "1.4.0",
#     keywords = ("pip", "single cell", "scRFE"),
#     description = "Single-cell identity definition using random forest modelling and recursive feature elimination",
#     long_description = "Single-cell identity one vs all classification using random forest modelling and recursive feature elimination",
#     license = "MIT Licence",
#
#     url = "https://github.com/czbiohub/scRFE",
#     author = "Madeline Park",
#     author_email = "madeline.park@czbiohub.org",
#
#     packages = find_packages(),
#     include_package_data = True,
#     platforms = "any",
#     install_requires=[
# 		'anndata>= 0.6.21',
#         'matplotlib>=3.1.1',
#         'numpy>=1.16.4',
#         'scikit-learn>=0.20.3',
#         'scanpy>=1.4.3',
#         'pandas>=0.24.2'
#
#     ]
# )
