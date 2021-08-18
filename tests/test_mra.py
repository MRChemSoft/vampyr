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

def test_BoundingBox_1D():
    world = vp.D1.BoundingBox(corner=[-1], nboxes=[2], scaling=[np.pi])
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

def test_BoundingBox_2D():
    world = vp.D2.BoundingBox(corner=[-1,-2], nboxes=[2, 4], scaling=[np.pi, np.pi])
    assert world.size() == 8
    assert world.size(dim=1) == 4
    assert world.getScale() == 0
    assert world.isPeriodic() == False
    assert world.getBoxLengths() == pytest.approx([2*np.pi, 4*np.pi])
    assert world.getBoxLength(dim=1) == pytest.approx(4*np.pi)
    assert world.getUpperBounds() == pytest.approx([np.pi, 2*np.pi])
    assert world.getUpperBound(dim=1) == pytest.approx(2*np.pi)
    assert world.getLowerBounds() == pytest.approx([-np.pi, -2*np.pi])
    assert world.getLowerBound(dim=1) == pytest.approx(-2*np.pi)
    assert world.getUnitLengths() == pytest.approx([np.pi, np.pi])
    assert world.getUnitLength(dim=1) == pytest.approx(np.pi)
    assert world.getScalingFactors() == pytest.approx([np.pi, np.pi])
    assert world.getScalingFactor(dim=1) == pytest.approx(np.pi)

def test_BoundingBox_3D():
    world = vp.D3.BoundingBox(corner=[-1,-2,-3], nboxes=[2, 4, 6], scaling=[np.pi, np.pi, np.pi])
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

def test_PeriodicBox_1D():
    world = vp.D1.BoundingBox(scaling=[np.pi], pbc=True)
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

def test_PeriodicBox_2D():
    world = vp.D2.BoundingBox(scaling=[np.pi, np.pi], pbc=True)
    assert world.size() == 1
    assert world.size(dim=1) == 1
    assert world.getScale() == 0
    assert world.isPeriodic() == True
    assert world.getBoxLengths() == pytest.approx([np.pi, np.pi])
    assert world.getBoxLength(dim=1) == pytest.approx(np.pi)
    assert world.getUpperBounds() == pytest.approx([np.pi, np.pi])
    assert world.getUpperBound(dim=1) == pytest.approx(np.pi)
    assert world.getLowerBounds() == pytest.approx([0.0, 0.0])
    assert world.getLowerBound(dim=1) == pytest.approx(0.0)
    assert world.getUnitLengths() == pytest.approx([np.pi, np.pi])
    assert world.getUnitLength(dim=1) == pytest.approx(np.pi)
    assert world.getScalingFactors() == pytest.approx([np.pi, np.pi])
    assert world.getScalingFactor(dim=1) == pytest.approx(np.pi)

def test_PeriodicBox_3D():
    world = vp.D3.BoundingBox(scaling=[np.pi, np.pi, np.pi], pbc=True)
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

def test_MRA_1D():
    legendre = vp.LegendreBasis(order=5)
    interpol = vp.InterpolatingBasis(order=5)
    world = vp.D1.BoundingBox(scale=-1)
    mra = vp.D1.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra == vp.D1.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra != vp.D1.MultiResolutionAnalysis(box=world, basis=legendre, max_depth=20)
    assert mra.getMaxDepth() == 20
    assert mra.getMaxScale() == 19
    assert mra.getWorldBox() == world
    assert mra.getScalingBasis() == interpol

def test_MRA_2D():
    legendre = vp.LegendreBasis(order=5)
    interpol = vp.InterpolatingBasis(order=5)
    world = vp.D2.BoundingBox(scale=-1)
    mra = vp.D2.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra == vp.D2.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra != vp.D2.MultiResolutionAnalysis(box=world, basis=legendre, max_depth=20)
    assert mra.getMaxDepth() == 20
    assert mra.getMaxScale() == 19
    assert mra.getWorldBox() == world
    assert mra.getScalingBasis() == interpol

def test_MRA_3D():
    legendre = vp.LegendreBasis(order=5)
    interpol = vp.InterpolatingBasis(order=5)
    world = vp.D3.BoundingBox(scale=-1)
    mra = vp.D3.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra == vp.D3.MultiResolutionAnalysis(box=world, basis=interpol, max_depth=20)
    assert mra != vp.D3.MultiResolutionAnalysis(box=world, basis=legendre, max_depth=20)
    assert mra.getMaxDepth() == 20
    assert mra.getMaxScale() == 19
    assert mra.getWorldBox() == world
    assert mra.getScalingBasis() == interpol
