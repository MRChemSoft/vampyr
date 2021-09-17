import numpy as np
from vampyr import vampyr3d as vp
import pytest

epsilon = 1.0e-3

D = 3
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
gauss = vp.GaussFunc(coef=alpha, exp=beta, pos=r0)

def test_Addition():
    tree_1 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    vp.project(prec=epsilon, out=tree_1, inp=gauss)

    tree_2 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_2, inp=tree_1)
    vp.add(out=tree_2, inp_a=tree_1, b=-0.5, inp_b=tree_1)
    assert tree_2.getNNodes() == tree_1.getNNodes()
    assert tree_2.integrate() == pytest.approx(0.5 * tree_1.integrate(), rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((-.5, tree_1))

    tree_3 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_3, inp=tree_vec_1)
    vp.add(out=tree_3, inp=tree_vec_1)
    assert tree_3.getNNodes() == tree_1.getNNodes()
    assert tree_3.integrate() == pytest.approx(0.5 * tree_1.integrate(), rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_4 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_4, inp=tree_vec_2)
    vp.add(out=tree_4, inp=tree_vec_2)
    assert tree_4.getNNodes() == tree_1.getNNodes()
    assert tree_4.integrate() == pytest.approx(2.0 * tree_1.integrate(), rel=epsilon)

    tree_5 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_5, inp=tree_1)
    tree_5.setZero()
    tree_5.add(inp=tree_1)
    tree_5.add(coef=1.0, inp=tree_1)
    assert tree_5.getNNodes() == tree_1.getNNodes()
    assert tree_5.integrate() == pytest.approx(2.0 * tree_1.integrate(), rel=epsilon)

def test_Multiplication():
    tree_1 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    vp.project(prec=epsilon, out=tree_1, inp=gauss)

    assert vp.dot(tree_1, tree_1) == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_2 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_2, inp=tree_1)
    vp.multiply(out=tree_2, inp_a=tree_1, inp_b=tree_1)
    assert tree_2.getNNodes() == tree_1.getNNodes()
    assert tree_2.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((1.0, tree_1))

    tree_3 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_3, inp=tree_vec_1)
    vp.multiply(out=tree_3, inp=tree_vec_1)
    assert tree_3.getNNodes() == tree_1.getNNodes()
    assert tree_3.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_4 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_4, inp=tree_vec_2)
    vp.multiply(out=tree_4, inp=tree_vec_2)
    assert tree_4.getNNodes() == tree_1.getNNodes()
    assert tree_4.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_5 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_5, inp=tree_1)
    vp.power(out=tree_5, inp=tree_1, pow=2.0)
    assert tree_5.getNNodes() == tree_1.getNNodes()
    assert tree_5.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_6 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_6, inp=tree_1)
    vp.square(out=tree_6, inp=tree_1)
    assert tree_6.getNNodes() == tree_1.getNNodes()
    assert tree_6.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_7 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_7, inp=tree_1)
    vp.copy_func(out=tree_7, inp=tree_1)
    tree_7.multiply(coef=1.0, inp=tree_1)
    assert tree_7.getNNodes() == tree_1.getNNodes()
    assert tree_7.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_8 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_8, inp=tree_1)
    vp.copy_func(out=tree_8, inp=tree_1)
    tree_8.square()
    assert tree_8.getNNodes() == tree_1.getNNodes()
    assert tree_8.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

    tree_9 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_9, inp=tree_1)
    vp.copy_func(out=tree_9, inp=tree_1)
    tree_9.power(pow=2.0)
    assert tree_9.getNNodes() == tree_1.getNNodes()
    assert tree_9.integrate() == pytest.approx(tree_1.getSquareNorm(), rel=epsilon)

def test_OverloadedOperators():
    tree_1 = vp.FunctionTree(mra)
    vp.build_grid(out=tree_1, inp=gauss)
    vp.project(out=tree_1, inp=gauss)
    ref_int = tree_1.integrate()
    ref_norm = tree_1.getSquareNorm()
    ref_nodes = tree_1.getNNodes()

    tree_2 = +tree_1
    assert tree_2.getNNodes() == ref_nodes
    assert tree_2.integrate() == pytest.approx(ref_int, rel=epsilon)

    tree_2 *= -2.0
    assert tree_2.getNNodes() == ref_nodes
    assert tree_2.integrate() == pytest.approx(-2.0*ref_int, rel=epsilon)

    tree_3 = -tree_1
    assert tree_3.getNNodes() == ref_nodes
    assert tree_3.integrate() == pytest.approx(-ref_int, rel=epsilon)

    tree_3 /= 2.0
    assert tree_3.getNNodes() == ref_nodes
    assert tree_3.integrate() == pytest.approx(-0.5*ref_int, rel=epsilon)

    tree_3 *= tree_2
    assert tree_3.getNNodes() > ref_nodes
    assert tree_3.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_4 = tree_1*2.0
    assert tree_4.getNNodes() == ref_nodes
    assert tree_4.integrate() == pytest.approx(2.0*ref_int, rel=epsilon)

    tree_5 = 2.0*tree_1
    assert tree_5.getNNodes() == ref_nodes
    assert tree_5.integrate() == pytest.approx(2.0*ref_int, rel=epsilon)

    tree_6 = tree_1/2.0 + 1.5*tree_1
    assert tree_6.getNNodes() == ref_nodes
    assert tree_6.integrate() == pytest.approx(2.0*ref_int, rel=epsilon)

    tree_6 += tree_1
    assert tree_6.getNNodes() == ref_nodes
    assert tree_6.integrate() == pytest.approx(3.0*ref_int, rel=epsilon)

    tree_6 -= 2.0*tree_4
    assert tree_6.getNNodes() == ref_nodes
    assert tree_6.integrate() == pytest.approx(-1.0*ref_int, rel=epsilon)

    tree_6 **= 2.0
    assert tree_6.getNNodes() > ref_nodes
    assert tree_6.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_7 = tree_1 - tree_2
    assert tree_7.getNNodes() == ref_nodes
    assert tree_7.integrate() == pytest.approx(3.0*ref_int, rel=epsilon)

    tree_8 = tree_1**2.0
    assert tree_8.getNNodes() > ref_nodes
    assert tree_8.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_9 = tree_1 * tree_1
    assert tree_9.getNNodes() > ref_nodes
    assert tree_9.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_10 = (tree_1 * tree_1).crop(epsilon) + (tree_1**2).crop(epsilon)
    assert tree_10.getNNodes() > tree_1.getNNodes()
    assert tree_10.getNNodes() < tree_9.getNNodes()
    assert tree_10.integrate() == pytest.approx(2.0*ref_norm, rel=epsilon)
