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


def version_from_file(rel_path):
    here = Path(__file__).parent
    f = here / rel_path
    with f.open("r") as fp:
        return fp.read().strip()


setup(
    version=version_from_file("VERSION"),
    packages=["vampyr"],
    package_dir={"": "src"},
    cmake_args=[
        "-DPYMOD_INSTALL_FULLDIR:STRING=src/vampyr",
    ],
)
