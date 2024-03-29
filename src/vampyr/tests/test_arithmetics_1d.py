import numpy as np
import pytest

from vampyr import vampyr1d as vp

epsilon = 1.0e-3

D = 1
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.8]
beta = 10.0
alpha = (beta / np.pi) ** (D / 2.0)
gauss = vp.GaussFunc(alpha=alpha, beta=beta, position=r0)


def test_Addition():
    tree_1 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_1, inp=gauss)
    vp.advanced.project(prec=epsilon, out=tree_1, inp=gauss)

    tree_2 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_2, inp=tree_1)
    vp.advanced.add(out=tree_2, inp_a=tree_1, b=-0.5, inp_b=tree_1)
    assert tree_2.nNodes() == tree_1.nNodes()
    assert tree_2.integrate() == pytest.approx(0.5 * tree_1.integrate(), rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((-0.5, tree_1))

    tree_3 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_3, inp=tree_vec_1)
    vp.advanced.add(out=tree_3, inp=tree_vec_1)
    assert tree_3.nNodes() == tree_1.nNodes()
    assert tree_3.integrate() == pytest.approx(0.5 * tree_1.integrate(), rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_4 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_4, inp=tree_vec_2)
    vp.advanced.add(out=tree_4, inp=tree_vec_2)
    assert tree_4.nNodes() == tree_1.nNodes()
    assert tree_4.integrate() == pytest.approx(2.0 * tree_1.integrate(), rel=epsilon)

    tree_5 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_5, inp=tree_1)
    tree_5.setZero()
    tree_5 += tree_1
    tree_5 += 1.0 * tree_1
    assert tree_5.nNodes() == tree_1.nNodes()
    assert tree_5.integrate() == pytest.approx(2.0 * tree_1.integrate(), rel=epsilon)


def test_Multiplication():
    tree_1 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_1, inp=gauss)
    vp.advanced.project(prec=epsilon, out=tree_1, inp=gauss)

    assert vp.dot(tree_1, tree_1) == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_2 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_2, inp=tree_1)
    vp.advanced.multiply(out=tree_2, inp_a=tree_1, inp_b=tree_1)
    assert tree_2.nNodes() == tree_1.nNodes()
    assert tree_2.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((1.0, tree_1))

    tree_3 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_3, inp=tree_vec_1)
    vp.advanced.multiply(out=tree_3, inp=tree_vec_1)
    assert tree_3.nNodes() == tree_1.nNodes()
    assert tree_3.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_4 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_4, inp=tree_vec_2)
    vp.advanced.multiply(out=tree_4, inp=tree_vec_2)
    assert tree_4.nNodes() == tree_1.nNodes()
    assert tree_4.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_5 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_5, inp=tree_1)
    vp.advanced.power(out=tree_5, inp=tree_1, pow=2.0)
    assert tree_5.nNodes() == tree_1.nNodes()
    assert tree_5.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_6 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_6, inp=tree_1)
    vp.advanced.square(out=tree_6, inp=tree_1)
    assert tree_6.nNodes() == tree_1.nNodes()
    assert tree_6.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_7 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_7, inp=tree_1)
    vp.advanced.copy_func(out=tree_7, inp=tree_1)
    tree_7 *= tree_1
    assert tree_7.nNodes() > tree_1.nNodes()
    assert tree_7.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_8 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_8, inp=tree_1)
    vp.advanced.copy_func(out=tree_8, inp=tree_1)
    tree_8 **= 2.0
    assert tree_8.nNodes() > tree_1.nNodes()
    assert tree_8.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)

    tree_9 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_9, inp=tree_1)
    vp.advanced.copy_func(out=tree_9, inp=tree_1)
    tree_9 **= 2.0
    assert tree_9.nNodes() > tree_1.nNodes()
    assert tree_9.integrate() == pytest.approx(tree_1.squaredNorm(), rel=epsilon)


def test_OverloadedOperators():
    tree_1 = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=tree_1, inp=gauss)
    vp.advanced.project(out=tree_1, inp=gauss)
    ref_int = tree_1.integrate()
    ref_norm = tree_1.squaredNorm()
    ref_nodes = tree_1.nNodes()

    tree_2 = +tree_1
    assert tree_2.nNodes() == ref_nodes
    assert tree_2.integrate() == pytest.approx(ref_int, rel=epsilon)

    tree_2 *= -2.0
    assert tree_2.nNodes() == ref_nodes
    assert tree_2.integrate() == pytest.approx(-2.0 * ref_int, rel=epsilon)

    tree_3 = -tree_1
    assert tree_3.nNodes() == ref_nodes
    assert tree_3.integrate() == pytest.approx(-ref_int, rel=epsilon)

    tree_3 /= 2.0
    assert tree_3.nNodes() == ref_nodes
    assert tree_3.integrate() == pytest.approx(-0.5 * ref_int, rel=epsilon)

    tree_3 *= tree_2
    assert tree_3.nNodes() > ref_nodes
    assert tree_3.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_4 = tree_1 * 2.0
    assert tree_4.nNodes() == ref_nodes
    assert tree_4.integrate() == pytest.approx(2.0 * ref_int, rel=epsilon)

    tree_5 = 2.0 * tree_1
    assert tree_5.nNodes() == ref_nodes
    assert tree_5.integrate() == pytest.approx(2.0 * ref_int, rel=epsilon)

    tree_6 = tree_1 / 2.0 + 1.5 * tree_1
    assert tree_6.nNodes() == ref_nodes
    assert tree_6.integrate() == pytest.approx(2.0 * ref_int, rel=epsilon)

    tree_6 += tree_1
    assert tree_6.nNodes() == ref_nodes
    assert tree_6.integrate() == pytest.approx(3.0 * ref_int, rel=epsilon)

    tree_6 -= 2.0 * tree_4
    assert tree_6.nNodes() == ref_nodes
    assert tree_6.integrate() == pytest.approx(-1.0 * ref_int, rel=epsilon)

    tree_6 **= 2.0
    assert tree_6.nNodes() > ref_nodes
    assert tree_6.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_7 = tree_1 - tree_2
    assert tree_7.nNodes() == ref_nodes
    assert tree_7.integrate() == pytest.approx(3.0 * ref_int, rel=epsilon)

    tree_8 = tree_1**2.0
    assert tree_8.nNodes() > ref_nodes
    assert tree_8.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_9 = tree_1 * tree_1
    assert tree_9.nNodes() > ref_nodes
    assert tree_9.integrate() == pytest.approx(ref_norm, rel=epsilon)

    tree_10 = (tree_1 * tree_1).crop(epsilon) + (tree_1**2).crop(epsilon)
    assert tree_10.nNodes() > tree_1.nNodes()
    assert tree_10.nNodes() < tree_9.nNodes()
    assert tree_10.integrate() == pytest.approx(2.0 * ref_norm, rel=epsilon)

    tree_vec_1 = []
    tree_vec_1.append((1.0, tree_1))
    tree_vec_1.append((-0.5, tree_1))

    tree_11 = vp.sum(tree_vec_1)
    assert tree_11.nNodes() == tree_1.nNodes()
    assert tree_11.integrate() == pytest.approx(0.5 * ref_int, rel=epsilon)

    tree_12 = vp.prod(tree_vec_1)
    assert tree_12.nNodes() > tree_1.nNodes()
    assert tree_12.integrate() == pytest.approx(-0.5 * ref_norm, rel=epsilon)

    tree_vec_2 = []
    tree_vec_2.append(tree_1)
    tree_vec_2.append(tree_1)

    tree_13 = vp.sum(tree_vec_2)
    assert tree_13.nNodes() == tree_1.nNodes()
    assert tree_13.integrate() == pytest.approx(2.0 * ref_int, rel=epsilon)

    tree_14 = vp.prod(tree_vec_2)
    assert tree_14.nNodes() > tree_1.nNodes()
    assert tree_14.integrate() == pytest.approx(ref_norm, rel=epsilon)
