import numpy as np
import vampyr as vp
import pytest

def test_LegendreBasis():
    basis = vp.LegendreBasis(order=5)
    assert basis != vp.InterpolatingBasis(order=5)
    assert basis != vp.LegendreBasis(order=6)
    assert basis == vp.LegendreBasis(order=5)
    assert basis.scalingType() == 0
    assert basis.scalingOrder() == 5
    assert basis.quadratureOrder() == 6

def test_InterpolatingBasis():
    basis = vp.InterpolatingBasis(order=5)
    assert basis != vp.LegendreBasis(order=5)
    assert basis != vp.InterpolatingBasis(order=6)
    assert basis == vp.InterpolatingBasis(order=5)
    assert basis.scalingType() == 1
    assert basis.scalingOrder() == 5
    assert basis.quadratureOrder() == 6

def test_InterpolatingScalingBasis():
    basis = vp.InterpolatingBasis(order=5)
    f = basis.scaling(i=2, l=0, n=1)
    dx = 0.0001
    x = np.arange(0, 0.51, dx)
    f_squared = [f([x])*f([x]) for x in x]
    assert 1.0 == pytest.approx(np.trapz(f_squared, dx=dx), 0.001)

def test_InterpolatingWaveletBasis():
    basis = vp.InterpolatingBasis(order=5)
    f = basis.wavelet(i=2, l=0, n=1)
    dx = 0.0001
    x = np.arange(0, 0.51, dx)
    f_squared = [f([x])*f([x]) for x in x]
    assert 1.0 == pytest.approx(np.trapz(f_squared, dx=dx), 0.001)

def test_LegendreScalingBasis():
    basis = vp.LegendreBasis(order=5)
    f = basis.scaling(i=2, l=0, n=1)
    dx = 0.0001
    x = np.arange(0, 0.51, dx)
    f_squared = [f([x])*f([x]) for x in x]
    assert 1.0 == pytest.approx(np.trapz(f_squared, dx=dx), 0.001)

def test_LegendreWaveletBasis():
    basis = vp.LegendreBasis(order=5)
    f = basis.wavelet(i=2, l=0, n=1)
    dx = 0.0001
    x = np.arange(0, 0.51, dx)
    f_squared = [f([x])*f([x]) for x in x]
    assert 1.0 == pytest.approx(np.trapz(f_squared, dx=dx), 0.001)
