![VAMPyR logo](https://github.com/MRChemSoft/VAMPyR/raw/master/doc/gfx/logo.png)

[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
[![Travis CI build status](https://travis-ci.org/MRChemSoft/vampyr.svg?branch=master)](https://travis-ci.org/MRChemSoft/vampyr)
[![CircleCI](https://circleci.com/gh/MRChemSoft/vampyr/tree/master.svg?style=svg)](https://circleci.com/gh/MRChemSoft/vampyr/tree/master)
[![codecov](https://codecov.io/gh/MRChemSoft/vampyr/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/vampyr)

The Very Accurate Multiresolution Python Routines (VAMPyR) package is a high
level Python interface to the [MRCPP](https://github.com/MRChemSoft/mrcpp) code.

### Installation

```sh
$ git clone https://github.com/MRChemSoft/vampyr.git
$ cd vampyr
$ git submodule update --init --recursive
$ ./setup
$ cd build
$ make
```

The vampyr.so file is now located in vampyr/build/lib, from there it can 
be copied into your python working directory.
