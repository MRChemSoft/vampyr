import numpy as np
import vampyr as vp
import pytest

epsilon = 1.0e-3

D = 3
s = 2
k = 5
N = -2
interpol = vp.InterpolatingBasis(order=k)
world = vp.D3.BoundingBox(scale=N)
mra = vp.D3.MultiResolutionAnalysis(box=world, basis=interpol)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
gauss = vp.D3.GaussFunc(coef=alpha, exp=beta, pos=r0)

def func(r):
    R2 = sum([(x1 - x2)**2 for x1, x2 in zip(r, r0)])
    return alpha*np.exp(-beta*R2)

def test_BuildProjectRefine():
    tree = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree, scales=s)
    assert tree.getDepth() == s + 1
    assert tree.getNEndNodes() == 2**(s*D)

    vp.project(out=tree, inp=func)
    assert tree.getDepth() == s + 1
    assert tree.getNEndNodes() == 2**(s*D)
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon);

    vp.refine_grid(out=tree, scales=1)
    assert tree.getDepth() == s + 2
    assert tree.getNEndNodes() == 2**((s+1)*D)
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon)

def test_BuildProjectCopyClear():
    tree_1 = vp.D3.FunctionTree(mra)
    tree_2 = vp.D3.FunctionTree(mra)

    vp.build_grid(out=tree_1, inp=gauss)
    assert tree_1.getDepth() > 1

    vp.project(prec=epsilon, out=tree_1, inp=gauss, abs_prec=True)
    assert tree_1.integrate() == pytest.approx(1.0, rel=epsilon);

    vp.copy_grid(out=tree_2, inp=tree_1)
    assert tree_2.getDepth() == tree_1.getDepth()

    vp.copy_func(out=tree_2, inp=tree_1)
    assert tree_2.integrate() == pytest.approx(tree_1.integrate(), rel=epsilon)
    assert tree_2.getSquareNorm() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    vp.clear_grid(out=tree_2)
    assert tree_2.getDepth() == tree_1.getDepth()
    assert tree_2.getSquareNorm() < 0.0

def test_BuildUnionGrid():
    # building sharp grid
    tree_1 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    assert tree_1.getDepth() > 1

    # building wide grid
    tree_2 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_2, scales=s)
    assert tree_2.getDepth() > 1

    # building union grid
    tree_3 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_3, inp=tree_1)
    vp.build_grid(out=tree_3, inp=tree_2)
    assert tree_3.getSquareNorm() < 0.0
    assert tree_3.getNNodes() > tree_1.getNNodes()
    assert tree_3.getNNodes() > tree_2.getNNodes()

    # refining union grid
    tree_4 = vp.D3.FunctionTree(mra)
    tree_4.setZero()
    while vp.refine_grid(out=tree_4, inp=tree_1) > 0:
        pass
    while vp.refine_grid(out=tree_4, inp=tree_2) > 0:
        pass
    assert tree_4.getSquareNorm() == pytest.approx(0.0, abs=epsilon)
    assert tree_4.getNNodes() == tree_3.getNNodes()

    # building union grid from vector
    tree_vec = []
    tree_vec.append(tree_1)
    tree_vec.append(tree_2)

    tree_5 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_5, inp=tree_vec)
    assert tree_5.getSquareNorm() < 0.0
    assert tree_5.getNNodes() == tree_3.getNNodes()

def test_ClearProjectRefine():
    tree = vp.D3.FunctionTree(mra)

    nsplit = 1
    while nsplit > 0:
        vp.clear_grid(out=tree)
        vp.project(out=tree, inp=func)
        nsplit = vp.refine_grid(out=tree, prec=epsilon)
        assert tree.getSquareNorm() > 0.0
    assert tree.integrate() == pytest.approx(1.0, rel=epsilon)
