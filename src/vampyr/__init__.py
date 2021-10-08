# -*- coding: utf-8 -*-

from ._vampyr import *

__version__ = _vampyr.__version__
__doc__ = _vampyr.__doc__

_dim_doc = """VAMPyR {dim:d}D
---------

{dim:d}-dimensional bindings to MRCPP.

.. currentmodule:: vampyr{dim:d}d

.. autosummary::
   :toctree: _generate
"""

vampyr1d.__doc__ = _dim_doc.format(dim=1)
vampyr2d.__doc__ = _dim_doc.format(dim=2)
vampyr3d.__doc__ = _dim_doc.format(dim=3)
