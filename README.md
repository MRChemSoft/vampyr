![VAMPyR logo](https://github.com/MRChemSoft/VAMPyR/raw/master/docs/gfx/logo.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4117602.svg)](https://doi.org/10.5281/zenodo.4117602)
[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
![Build and test VAMPyR](https://github.com/MRChemSoft/vampyr/workflows/Build%20and%20test%20VAMPyR/badge.svg)
[![codecov](https://codecov.io/gh/MRChemSoft/vampyr/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/vampyr)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MRChemSoft/vampyr/master?urlpath=lab%2Ftree%2Fdocs%2Fnotebooks)

The Very Accurate Multiresolution Python Routines (VAMPyR) package is a high
level Python interface to the [MRCPP](https://github.com/MRChemSoft/mrcpp) code.

## Installation

### From source

To build VAMPyR from source with OpenMP parallelization using CMake:

```sh
$ git clone https://github.com/MRChemSoft/vampyr.git
$ cd vampyr
$ mkdir build
$ cd build
$ cmake ..
$ make
```

To use the VAMPyR module you need to add it to your `$PYTHONPATH`:

``` sh
$ export PYTHONPATH=<path-to-vampyr>/build/lib/python<version>/site-packages:$PYTHONPATH
```

### Using pip

To build VAMPyR from source with OpenMP parallelization using `pip`, here in a virtual
environment called `myenv`:

```sh
$ git clone https://github.com/MRChemSoft/vampyr.git
$ cd vampyr
$ virtualenv myenv
$ source myenv/bin/activate
$ python -m pip install .
```

The VAMPyR module is now available whenever you have activated the `myenv` environment.

### Using Conda

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/vampyr/badges/version.svg)](https://anaconda.org/conda-forge/vampyr)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/vampyr/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/vampyr)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/vampyr/badges/downloads.svg)](https://anaconda.org/conda-forge/vampyr)

To install VAMPyR in a Conda environment `myenv`:

    $ conda create -n myenv
    $ conda activate myenv
    $ conda install -c conda-forge vampyr               # latest version (OpenMP)
    $ conda install -c conda-forge vampyr=0.1.0rc0      # tagged version (OpenMP)

To list all available versions:

    $ conda search -c conda-forge vampyr

Note that the conda-forge package is _always_ built with OpenMP support enabled
in the MRCPP backend.

The VAMPyR module is now available whenever you have activated the `myenv` environment.

### Creating a Conda environment from a .yml file

You can also create a Conda environment from a .yml file that already specifies VAMPyR and
other useful packages such as numpy, and matplotlib. Here's how:

1. Write an `environment.yml` file, for example:

    ```yaml
    name: myenv
    channels:
      - conda-forge
    dependencies:
      - vampyr
      - numpy
      - matplotlib
      - jupyterlab
    ```

2. Create the environment from the `environment.yml` file:

    ```sh
    $ conda env create -f environment.yml
    ```
3. Activate the environment:

    ```sh
    $ conda activate myenv
    ```
The VAMPyR module, along with numpy and matplotlib, is now available whenever
you have activated the myenv environment.