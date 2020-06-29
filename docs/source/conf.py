# .readthedocs.yml
# Read the Docs configuration file
# See https://docs.readthedocs.io/en/stable/config-file/v2.html for details

# Required
version: 2

# Build documentation in the docs/ directory with Sphinx
sphinx:
  configuration: docs/conf.py

# Build documentation with MkDocs
#mkdocs:
#  configuration: mkdocs.yml

# Optionally build your docs in additional formats such as PDF
formats:
  - pdf

# Optionally set the version of Python and requirements required to build your docs
# python:
#   version: 3.7
#   install:
#     - requirements: docs/requirements.txt
