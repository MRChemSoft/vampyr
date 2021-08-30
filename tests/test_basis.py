import numpy as np
import vampyr as vp
import pytest

def test_LegendreBasis():
    basis = vp.LegendreBasis(order=5)
    assert basis != vp.InterpolatingBasis(order=5)
    assert basis != vp.LegendreBasis(order=6)
    assert basis == vp.LegendreBasis(order=5)
    assert basis.getScalingType() == 0
    assert basis.getScalingOrder() == 5
    assert basis.getQuadratureOrder() == 6

def test_InterpolatingBasis():
    basis = vp.InterpolatingBasis(order=5)
    assert basis != vp.LegendreBasis(order=5)
    assert basis != vp.InterpolatingBasis(order=6)
    assert basis == vp.InterpolatingBasis(order=5)
    assert basis.getScalingType() == 1
    assert basis.getScalingOrder() == 5
    assert basis.getQuadratureOrder() == 6
