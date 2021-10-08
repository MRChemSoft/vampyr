# -*- coding: utf-8 -*-

import sys
from pathlib import Path

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages


def version_from_file(rel_path):
    here = Path(__file__).parent
    f = here / rel_path
    with f.open("r") as fp:
        return fp.read().strip()


setup(
    name="vampyr",
    version=version_from_file("VERSION"),
    description="The Very Accurate Multiresolution Python Routines (VAMPyR) package is a high level Python interface to the MRCPP code.",
    author="Magnar Bjørgve",
    license="LGPL",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_args=[
        "-DCMAKE_JOB_POOL_COMPILE:STRING=compile",
        "-DCMAKE_JOB_POOL_LINK:STRING=link",
        "-DCMAKE_JOB_POOLS:STRING=compile=4;link=2",
    ],
    include_package_data=True,
    extras_require={"test": ["numpy", "pytest"]},
)
