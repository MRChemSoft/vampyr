import pytest

from vampyr import vampyr3d as vp

def test_ScalingProjector():
    def f(x):
        return x

    mra = vp.MultiResolutionAnalysis(box=[0, 1], order=7)
    P_scaling = vp.ScalingProjector(mra, 2)
    P_wavelet = vp.WaveletProjector(mra, 2)

    with pytest.raises(Exception):
        P_scaling(f)

    with pytest.raises(Exception):
        P_wavelet(f)