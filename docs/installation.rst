.. _installation:
============
Installation
============


Prerequisites
=============
We recommend to use an environment manager like conda to install **VAMPyR** and its prerequisites.

To set up an environment named `vamp-env` and activating it using conda we can run:
```
conda create -n vamp-env
conda activate vamp-env
```
Other environment managers are also possible, but we will not cover them here.


We will work from this environment for the rest of the installation process.

Install options
===============

From conda
----------

The easiest way to install **VAMPyR** is to use the conda package manager.
You can install **VAMPyR** by running the following commmand while in the `vamp-env` environment:
```bash
conda install -c conda-forge vampyr
```
This should install **VAMPyR** and all its dependencies. To use **VAMPyR** in a script or a Jupyter notebook, you can simply follow the instructions in the notebooks.


From source
-----------

You can get the source code from github by running:
```bash
git clone https://github.com/MRChemSoft/vampyr.git
```
This will create a directory called `vampyr` in your current working directory, and download the source code there.

Alternatively one can run
```bash
git clone git@github.com:MRChemSoft/vampyr.git
```
to clone the repository using ssh. This requires that you have set up an ssh key with github.



Here we need to write about how to install **VAMPyR**. Possibilities:
conda package, build from source, ....


