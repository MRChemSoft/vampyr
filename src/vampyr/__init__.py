# -*- coding: utf-8 -*-

from ._vampyr import *
from .environ import _set_mwfilters_path

__version__ = _vampyr.__version__
__doc__ = _vampyr.__doc__

_dim_doc = """VAMPyR {dim:d}D
---------

{dim:d}-dimensional bindings to MRCPP.

.. currentmodule:: vampyr{dim:d}d

.. autosummary::
   :toctree: generate
"""

vampyr1d.__doc__ = _dim_doc.format(dim=1)
vampyr2d.__doc__ = _dim_doc.format(dim=2)
vampyr3d.__doc__ = _dim_doc.format(dim=3)

_REAL_DTYPES = (float, "float", "double", None)
_COMPLEX_DTYPES = (complex, "complex")


def _install_dtype_projector(sub):
    """Wrap the native projectors in a dtype-dispatching factory.

    ``ScalingProjector(mra, prec, dtype=complex)`` returns a projector that
    builds ``ComplexFunctionTree``s from complex-valued callables, while the
    default ``dtype=float`` preserves the original real behaviour, including
    the ``prec``/``scale`` constructor overloads.
    """
    _real_cls = sub.ScalingProjector
    _complex_cls = sub.ComplexScalingProjector

    def ScalingProjector(mra, *args, dtype=float, **kwargs):
        if dtype in _REAL_DTYPES:
            cls = _real_cls
        elif dtype in _COMPLEX_DTYPES:
            cls = _complex_cls
        else:
            raise TypeError(f"Unsupported dtype for ScalingProjector: {dtype!r}")
        return cls(mra, *args, **kwargs)

    ScalingProjector.__doc__ = (
        "Scaling projector onto the MRA. Accepts either prec (adaptive) or "
        "scale (fixed), and dtype=float (default) or dtype=complex."
    )
    sub._RealScalingProjector = _real_cls
    sub.ScalingProjector = ScalingProjector


def _install_dtype_function_map(sub):
    """Wrap the native function maps in a dtype-dispatching factory.

    ``FunctionMap(fmap, prec, dtype=complex)`` returns a map that consumes and
    produces ``ComplexFunctionTree``s, while the default ``dtype=float``
    preserves the original real behaviour. Mirrors
    ``_install_dtype_projector`` above; both the real and the complex class
    take the same ``(fmap, prec)`` constructor, so no ``*args`` juggling is
    needed here.
    """
    _real_cls = sub.FunctionMap
    _complex_cls = sub.ComplexFunctionMap

    def FunctionMap(fmap, prec, dtype=float, **kwargs):
        if dtype in _REAL_DTYPES:
            cls = _real_cls
        elif dtype in _COMPLEX_DTYPES:
            cls = _complex_cls
        else:
            raise TypeError(f"Unsupported dtype for FunctionMap: {dtype!r}")
        return cls(fmap, prec, **kwargs)

    FunctionMap.__doc__ = (
        "Pointwise map applied through the MW representation. Takes the map "
        "and a precision, and dtype=float (default) or dtype=complex."
    )
    sub._RealFunctionMap = _real_cls
    sub.FunctionMap = FunctionMap


for _sub in (vampyr1d, vampyr2d, vampyr3d):
    _install_dtype_projector(_sub)
    _install_dtype_function_map(_sub)
del _sub

_set_mwfilters_path()
