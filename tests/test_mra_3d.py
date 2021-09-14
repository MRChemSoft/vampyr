import numpy as np
from vampyr import vampyr3d as vp
from vampyr import InterpolatingBasis, LegendreBasis
import pytest

def test_BoundingBox():
    world = vp.BoundingBox(corner=[-1,-2,-3], nboxes=[2, 4, 6], scaling=[np.pi, np.pi, np.pi])
    assert world.size() == 48
    assert world.size(dim=2) == 6
    assert world.getScale() == 0
    assert world.isPeriodic() == False
    assert world.getBoxLengths() == pytest.approx([2*np.pi, 4*np.pi, 6*np.pi])
    assert world.getBoxLength(dim=2) == pytest.approx(6*np.pi)
    assert world.getUpperBounds() == pytest.approx([np.pi, 2*np.pi, 3*np.pi])
    assert world.getUpperBound(dim=2) == pytest.approx(3*np.pi)
    assert world.getLowerBounds() == pytest.approx([-np.pi, -2*np.pi, -3*np.pi])
    assert world.getLowerBound(dim=2) == pytest.approx(-3*np.pi)
    assert world.getUnitLengths() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.getUnitLength(dim=2) == pytest.approx(np.pi)
    assert world.getScalingFactors() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.getScalingFactor(dim=2) == pytest.approx(np.pi)

def test_PeriodicBox():
    world = vp.BoundingBox(scaling=[np.pi, np.pi, np.pi], pbc=True)
    assert world.size() == 1
    assert world.size(dim=1) == 1
    assert world.getScale() == 0
    assert world.isPeriodic() == True
    assert world.getBoxLengths() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.getBoxLength(dim=1) == pytest.approx(np.pi)
    assert world.getUpperBounds() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.getUpperBound(dim=1) == pytest.approx(np.pi)
    assert world.getLowerBounds() == pytest.approx([0.0, 0.0, 0.0])
    assert world.getLowerBound(dim=1) == pytest.approx(0.0)
    assert world.getUnitLengths() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.getUnitLength(dim=1) == pytest.approx(np.pi)
    assert world.getScalingFactors() == pytest.approx([np.pi, np.pi, np.pi])
    assert world.getScalingFactor(dim=1) == pytest.approx(np.pi)

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
