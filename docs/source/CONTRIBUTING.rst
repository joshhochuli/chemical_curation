.. highlight:: shell

============
CONTRIBUTING
============

Getting Setup
------------

This section will get you set up with the conda build process.

1. Fork the `chemical_curation` repo on GitHub.
2. Clone your fork locally::

     $ git clone git@github.com:your_github_username/chemical_curation.git && cd chemical_curation

3. Install `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>` if you
   do not already have it.
4. Create a new conda environment and install conda-build and conda-verify::

     $ conda create --name tbv_dev conda-build conda-verify
     $ conda activate tbv_dev

5. The project can now be built locally using::

     $ conda-build -c conda-forge .

6. Near the end of the conda-build output, there should a line that says "TEST
   END: " with the location of the package. Copy the path up to the last
   directory name (so if the full output is TEST END: /path/to/pkg/pkg.tar.bz2,
   copy /big/to/pkg/). Use this directory as the channel for conda install::

     $ conda install -c /path/to/pkg/ chemical_curation

   You should now be able to import modules from the package, as long as you are
   in the conda environment where you built it.
     
