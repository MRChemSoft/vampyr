import numpy as np
from vampyr import vampyr1d as vp
from vampyr import InterpolatingBasis, LegendreBasis
import pytest

def test_BoundingBox():
    world = vp.BoundingBox(corner=[-1], nboxes=[2], scaling=[np.pi])
    assert world.size() == 2
    assert world.size(dim=0) == 2
    assert world.getScale() == 0
    assert world.isPeriodic() == False
    assert world.getBoxLengths() == pytest.approx([2*np.pi])
    assert world.getBoxLength(dim=0) == pytest.approx(2*np.pi)
    assert world.getUpperBounds() == pytest.approx([np.pi])
    assert world.getUpperBound(dim=0) == pytest.approx(np.pi)
    assert world.getLowerBounds() == pytest.approx([-np.pi])
    assert world.getLowerBound(dim=0) == pytest.approx(-np.pi)
    assert world.getUnitLengths() == pytest.approx([np.pi])
    assert world.getUnitLength(dim=0) == pytest.approx(np.pi)
    assert world.getScalingFactors() == pytest.approx([np.pi])
    assert world.getScalingFactor(dim=0) == pytest.approx(np.pi)

def test_PeriodicBox():
    world = vp.BoundingBox(scaling=[np.pi], pbc=True)
    assert world.size() == 1
    assert world.size(dim=0) == 1
    assert world.getScale() == 0
    assert world.isPeriodic() == True
    assert world.getBoxLengths() == pytest.approx([np.pi])
    assert world.getBoxLength(dim=0) == pytest.approx(np.pi)
    assert world.getUpperBounds() == pytest.approx([np.pi])
    assert world.getUpperBound(dim=0) == pytest.approx(np.pi)
    assert world.getLowerBounds() == pytest.approx([0.0])
    assert world.getLowerBound(dim=0) == pytest.approx(0.0)
    assert world.getUnitLengths() == pytest.approx([np.pi])
    assert world.getUnitLength(dim=0) == pytest.approx(np.pi)
    assert world.getScalingFactors() == pytest.approx([np.pi])
    assert world.getScalingFactor(dim=0) == pytest.approx(np.pi)

def test_MRA():
    legendre = LegendreBasis(order=5)
    interpol = InterpolatingBasis(order=5)
    world = vp.BoundingBox(scale=-1)
    mra = vp.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra == vp.MultiResolutionAnalysis(box=world, order=5, max_depth=20)
    assert mra == vp.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra != vp.MultiResolutionAnalysis(box=world, basis=legendre, max_depth=20)
    assert mra.getMaxDepth() == 20
    assert mra.getMaxScale() == 19
    assert mra.getWorldBox() == world
    assert mra.getScalingBasis() == interpol
