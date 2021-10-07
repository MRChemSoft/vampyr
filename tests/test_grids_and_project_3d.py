import numpy as np
from vampyr import vampyr3d as vp
import pytest

epsilon = 1.0e-3

D = 3
s = 2
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
gauss = vp.GaussFunc(coef=alpha, exp=beta, pos=r0)

def func(r):
    R2 = sum([(x1 - x2)**2 for x1, x2 in zip(r, r0)])
    return alpha*np.exp(-beta*R2)

def test_BuildProjectRefineCrop():
    tree = vp.FunctionTree(mra)
    vp.build_grid(out=tree, scales=s)
    assert tree.depth() == s + 1
    assert tree.nEndNodes() == 2**(s*D)

    vp.project(out=tree, inp=func)
    assert tree.depth() == s + 1
    assert tree.nEndNodes() == 2**(s*D)
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon);

    vp.refine_grid(out=tree, scales=1)
    assert tree.depth() == s + 2
    assert tree.nEndNodes() == 2**((s+1)*D)
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon)

    pre_crop = tree.nNodes()
    tree.crop(prec=epsilon)
    post_crop = tree.nNodes()
    assert pre_crop > post_crop

def test_BuildProjectCopyClear():
    tree_1 = vp.FunctionTree(mra)
    tree_2 = vp.FunctionTree(mra)

    vp.build_grid(out=tree_1, inp=gauss)
    assert tree_1.depth() > 1

    vp.project(prec=epsilon, out=tree_1, inp=gauss, abs_prec=True)
    assert tree_1.integrate() == pytest.approx(1.0, rel=epsilon);

    vp.copy_grid(out=tree_2, inp=tree_1)
    assert tree_2.depth() == tree_1.depth()

    vp.copy_func(out=tree_2, inp=tree_1)
    assert tree_2.integrate() == pytest.approx(tree_1.integrate(), rel=epsilon)
    assert tree_2.norm() == pytest.approx(tree_1.norm(), rel=epsilon)

    vp.clear_grid(out=tree_2)
    assert tree_2.depth() == tree_1.depth()
    assert tree_2.norm() < 0.0

def test_BuildUnionGrid():
    # building sharp grid
    tree_1 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    assert tree_1.depth() > 1

    # building wide grid
    tree_2 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_2, scales=s)
    assert tree_2.depth() > 1

    # building union grid
    tree_3 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_3, inp=tree_1)
    vp.build_grid(out=tree_3, inp=tree_2)
    assert tree_3.norm() < 0.0
    assert tree_3.nNodes() > tree_1.nNodes()
    assert tree_3.nNodes() > tree_2.nNodes()

    # refining union grid
    tree_4 = vp.FunctionTree(mra)
    tree_4.setZero()
    while vp.refine_grid(out=tree_4, inp=tree_1) > 0:
        pass
    while vp.refine_grid(out=tree_4, inp=tree_2) > 0:
        pass
    assert tree_4.norm() == pytest.approx(0.0, abs=epsilon)
    assert tree_4.nNodes() == tree_3.nNodes()

    # building union grid from vector
    tree_vec = []
    tree_vec.append(tree_1)
    tree_vec.append(tree_2)

    tree_5 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_5, inp=tree_vec)
    assert tree_5.norm() < 0.0
    assert tree_5.nNodes() == tree_3.nNodes()

def test_ClearProjectRefine():
    tree = vp.FunctionTree(mra)

    nsplit = 1
    while nsplit > 0:
        vp.clear_grid(out=tree)
        vp.project(out=tree, inp=func)
        nsplit = vp.refine_grid(out=tree, prec=epsilon)
        assert tree.norm() > 0.0
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon)

def test_ProjectRescaleCopyNormalize():
    Q = vp.MWProjector(mra, epsilon)
    tree = Q(gauss)
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon)

    tree.rescale(coef=np.pi)
    assert tree.integrate() == pytest.approx(np.pi, rel=epsilon)

    tree_shallow = tree
    tree_deep = tree.deepCopy()
    assert tree_shallow.integrate() == pytest.approx(np.pi, rel=epsilon)
    assert tree_deep.integrate() == pytest.approx(np.pi, rel=epsilon)

    tree.normalize()
    assert tree.norm() == pytest.approx(1.0, rel=epsilon)
    assert tree_deep.integrate() == pytest.approx(np.pi, rel=epsilon)
    assert tree_shallow.integrate() != pytest.approx(np.pi, rel=epsilon)

def test_BuildProjectSemiPeriodicGauss():
    sfac = [np.pi/3, np.pi/3, np.pi/3]
    periodic_world = vp.BoundingBox(scaling=sfac, pbc=True)
    pbc = vp.MultiResolutionAnalysis(box=periodic_world, order=k)

    tree_1 = vp.FunctionTree(pbc)
    vp.build_grid(out=tree_1, inp=gauss)
    vp.project(out=tree_1, inp=gauss)
    assert tree_1.integrate() < 1.0

    pgauss = gauss.periodify(period=sfac, std_dev=6.0)
    assert pgauss.size() > 0

    tree_2 = vp.FunctionTree(pbc)
    vp.build_grid(out=tree_2, inp=pgauss)
    vp.project(out=tree_2, inp=pgauss)
    assert tree_2.integrate() == pytest.approx(1.0, rel=epsilon)

def test_BuildProjectSemiPeriodicGaussExp():
    gexp = vp.GaussExp()
    gexp.append(gauss)
    gexp.append(gauss)

    sfac = [np.pi/3, np.pi/3, np.pi/3]
    periodic_world = vp.BoundingBox(scaling=sfac, pbc=True)
    pbc = vp.MultiResolutionAnalysis(box=periodic_world, order=k)

    tree_1 = vp.FunctionTree(pbc)
    vp.build_grid(out=tree_1, inp=gexp)
    vp.project(out=tree_1, inp=gexp)
    assert tree_1.integrate() < 2.0

    pexp = gexp.periodify(period=sfac, std_dev=6.0)
    assert pexp.size() > 2

    tree_2 = vp.FunctionTree(pbc)
    vp.build_grid(out=tree_2, inp=pexp)
    vp.project(out=tree_2, inp=pexp)
    assert tree_2.integrate() == pytest.approx(2.0, rel=epsilon)
