import numpy as np
import pytest

from vampyr import InterpolatingBasis, LegendreBasis
from vampyr import vampyr1d as vp


def test_BoundingBox():
    world = vp.BoundingBox(corner=[-1], nboxes=[2], scaling=[np.pi])
    assert world.size() == 2
    assert world.size(dim=0) == 2
    assert world.scale() == 0
    assert world.isPeriodic() == False
    assert world.boxLengths() == pytest.approx([2 * np.pi])
    assert world.boxLength(dim=0) == pytest.approx(2 * np.pi)
    assert world.upperBounds() == pytest.approx([np.pi])
    assert world.upperBound(dim=0) == pytest.approx(np.pi)
    assert world.lowerBounds() == pytest.approx([-np.pi])
    assert world.lowerBound(dim=0) == pytest.approx(-np.pi)
    assert world.unitLengths() == pytest.approx([np.pi])
    assert world.unitLength(dim=0) == pytest.approx(np.pi)
    assert world.scalingFactors() == pytest.approx([np.pi])
    assert world.scalingFactor(dim=0) == pytest.approx(np.pi)


def test_PeriodicBox():
    world = vp.BoundingBox(scaling=[np.pi], pbc=True)
    assert world.size() == 1
    assert world.size(dim=0) == 1
    assert world.scale() == 0
    assert world.isPeriodic() == True
    assert world.boxLengths() == pytest.approx([np.pi])
    assert world.boxLength(dim=0) == pytest.approx(np.pi)
    assert world.upperBounds() == pytest.approx([np.pi])
    assert world.upperBound(dim=0) == pytest.approx(np.pi)
    assert world.lowerBounds() == pytest.approx([0.0])
    assert world.lowerBound(dim=0) == pytest.approx(0.0)
    assert world.unitLengths() == pytest.approx([np.pi])
    assert world.unitLength(dim=0) == pytest.approx(np.pi)
    assert world.scalingFactors() == pytest.approx([np.pi])
    assert world.scalingFactor(dim=0) == pytest.approx(np.pi)


def test_MRA():
    legendre = LegendreBasis(order=5)
    interpol = InterpolatingBasis(order=5)
    world = vp.BoundingBox(scale=-1)
    mra = vp.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra == vp.MultiResolutionAnalysis(box=[0,2], order=5, max_depth=20)
    assert mra == vp.MultiResolutionAnalysis(box=world, order=5, max_depth=20)
    assert mra != vp.MultiResolutionAnalysis(box=world, basis=legendre, max_depth=20)
    assert mra.maxDepth() == 20
    assert mra.maxScale() == 19
    assert mra.world() == world
    assert mra.basis() == interpol
