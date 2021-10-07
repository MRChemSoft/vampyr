import numpy as np
from vampyr import vampyr3d as vp
from vampyr import InterpolatingBasis, LegendreBasis
import pytest

def test_BoundingBox():
    world = vp.BoundingBox(corner=[-1,-2,-3], nboxes=[2, 4, 6], scaling=[np.pi, np.pi, np.pi])
    assert world.size() == 48
    assert world.size(dim=2) == 6
    assert world.scale() == 0
    assert world.isPeriodic() == False
    assert world.boxLengths() == pytest.approx([2*np.pi, 4*np.pi, 6*np.pi])
    assert world.boxLength(dim=2) == pytest.approx(6*np.pi)
    assert world.upperBounds() == pytest.approx([np.pi, 2*np.pi, 3*np.pi])
    assert world.upperBound(dim=2) == pytest.approx(3*np.pi)
    assert world.lowerBounds() == pytest.approx([-np.pi, -2*np.pi, -3*np.pi])
    assert world.lowerBound(dim=2) == pytest.approx(-3*np.pi)
    assert world.unitLengths() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.unitLength(dim=2) == pytest.approx(np.pi)
    assert world.scalingFactors() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.scalingFactor(dim=2) == pytest.approx(np.pi)

def test_PeriodicBox():
    world = vp.BoundingBox(scaling=[np.pi, np.pi, np.pi], pbc=True)
    assert world.size() == 1
    assert world.size(dim=1) == 1
    assert world.scale() == 0
    assert world.isPeriodic() == True
    assert world.boxLengths() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.boxLength(dim=1) == pytest.approx(np.pi)
    assert world.upperBounds() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.upperBound(dim=1) == pytest.approx(np.pi)
    assert world.lowerBounds() == pytest.approx([0.0, 0.0, 0.0])
    assert world.lowerBound(dim=1) == pytest.approx(0.0)
    assert world.unitLengths() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.unitLength(dim=1) == pytest.approx(np.pi)
    assert world.scalingFactors() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.scalingFactor(dim=1) == pytest.approx(np.pi)

def test_MRA():
    legendre = LegendreBasis(order=5)
    interpol = InterpolatingBasis(order=5)
    world = vp.BoundingBox(scale=-1)
    mra = vp.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra == vp.MultiResolutionAnalysis(box=world, order=5, max_depth=20)
    assert mra == vp.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra != vp.MultiResolutionAnalysis(box=world, basis=legendre, max_depth=20)
    assert mra.maxDepth() == 20
    assert mra.maxScale() == 19
    assert mra.world() == world
    assert mra.basis() == interpol
