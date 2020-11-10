![VAMPyR logo](https://github.com/MRChemSoft/VAMPyR/raw/master/doc/gfx/logo.png)

[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
![Build and test VAMPyR](https://github.com/MRChemSoft/vampyr/workflows/Build%20and%20test%20VAMPyR/badge.svg)
[![codecov](https://codecov.io/gh/MRChemSoft/vampyr/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/vampyr)

The Very Accurate Multiresolution Python Routines (VAMPyR) package is a high
level Python interface to the [MRCPP](https://github.com/MRChemSoft/mrcpp) code.

### Installation

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
