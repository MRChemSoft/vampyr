import numpy as np
import vampyr as vp
import pytest

epsilon = 1.0e-3

D = 3
k = 5
N = -2
interpol = vp.InterpolatingBasis(order=k)
world = vp.D3.BoundingBox(scale=N)
mra = vp.D3.MultiResolutionAnalysis(box=world, basis=interpol)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
gauss = vp.D3.GaussFunc(coef=alpha, exp=beta, pos=r0)

def test_Addition():
    tree_1 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    vp.project(prec=epsilon, out=tree_1, inp=gauss)

    tree_2 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_2, inp=tree_1)
    vp.add(out=tree_2, inp_a=tree_1, b=-0.5, inp_b=tree_1)
    assert tree_2.getNNodes() == tree_1.getNNodes()
    assert tree_2.integrate() == pytest.approx(0.5 * tree_1.integrate(), rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((-.5, tree_1))

    tree_3 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_3, inp=tree_vec_1)
    vp.add(out=tree_3, inp=tree_vec_1)
    assert tree_3.getNNodes() == tree_1.getNNodes()
    assert tree_3.integrate() == pytest.approx(0.5 * tree_1.integrate(), rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_4 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_4, inp=tree_vec_2)
    vp.add(out=tree_4, inp=tree_vec_2)
    assert tree_4.getNNodes() == tree_1.getNNodes()
    assert tree_4.integrate() == pytest.approx(2.0 * tree_1.integrate(), rel=epsilon)

def test_Multiplication():
    tree_1 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    vp.project(prec=epsilon, out=tree_1, inp=gauss)

    assert vp.dot(tree_1, tree_1) == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_2 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_2, inp=tree_1)
    vp.multiply(out=tree_2, inp_a=tree_1, inp_b=tree_1)
    assert tree_2.getNNodes() == tree_1.getNNodes()
    assert tree_2.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((1.0, tree_1))

    tree_3 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_3, inp=tree_vec_1)
    vp.multiply(out=tree_3, inp=tree_vec_1)
    assert tree_3.getNNodes() == tree_1.getNNodes()
    assert tree_3.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_4 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_4, inp=tree_vec_2)
    vp.multiply(out=tree_4, inp=tree_vec_2)
    assert tree_4.getNNodes() == tree_1.getNNodes()
    assert tree_4.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_5 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_5, inp=tree_1)
    vp.power(out=tree_5, inp=tree_1, pow=2.0)
    assert tree_5.getNNodes() == tree_1.getNNodes()
    assert tree_5.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_6 = vp.D3.FunctionTree(mra)
    vp.build_grid(out=tree_6, inp=tree_1)
    vp.square(out=tree_6, inp=tree_1)
    assert tree_6.getNNodes() == tree_1.getNNodes()
    assert tree_6.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)
