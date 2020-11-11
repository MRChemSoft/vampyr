![VAMPyR logo](https://github.com/MRChemSoft/VAMPyR/raw/master/doc/gfx/logo.png)

[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
![Build and test VAMPyR](https://github.com/MRChemSoft/vampyr/workflows/Build%20and%20test%20VAMPyR/badge.svg)
[![Travis CI build status](https://travis-ci.org/MRChemSoft/vampyr.svg?branch=master)](https://travis-ci.org/MRChemSoft/vampyr)
[![codecov](https://codecov.io/gh/MRChemSoft/vampyr/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/vampyr)

The Very Accurate Multiresolution Python Routines (VAMPyR) package is a high
level Python interface to the [MRCPP](https://github.com/MRChemSoft/mrcpp) code.

## Installation

### From source

To build VAMPyR from source with OpenMP parallelization:

```sh
$ git clone https://github.com/MRChemSoft/vampyr.git
$ cd vampyr
$ ./setup
$ cd build
$ make
```

To use the Python module:

``` sh
$ export PYTHONPATH=build/lib/python:$PYTHONPATH
```

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
