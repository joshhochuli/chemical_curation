.. highlight:: shell

============
CONTRIBUTING
============

Getting Setup
------------

This section will get you set up with the code in devlopment mode so you can
easily make and test changes.

1. Fork the `chemical_curation` repo on GitHub.
2. Clone your fork locally::

     $ git clone git@github.com:your_github_username/chemical_curation.git && cd chemical_curation

3. Install `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>` if you
   do not already have it.
4. Create an environment using the `environment.yaml` file in this repo. This
   will create a virtual environment named "chemical_curation"::

     $ conda env create -f environment.yaml
     $ conda activate chemical_curation

5. Run `pip install --editable .` to install the package in development mode;
   this means that changes you make to the code will be automatically applied to
   your local version of the package so that you don't have to reinstall after
   every change.   
