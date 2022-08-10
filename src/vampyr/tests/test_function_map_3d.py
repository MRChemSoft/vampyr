import numpy as np
import pytest

from vampyr import vampyr3d as vp

L = 1
k = 5
epsilon = 1.0e-3
mra = vp.MultiResolutionAnalysis(box=[-L, L], order=k)


def f(x):
    return x[0]


def g(x):
    return x[0] ** 2


def gmap(x):
    return x**2


P = vp.ScalingProjector(mra, prec=epsilon)
f_tree = P(f)
g_ref = P(g)

r_0 = [0.55, 0.45, 0.35]


def test_FunctionMap():
    g_tree = vp.FunctionTree(mra)
    vp.advanced.map(prec=epsilon, out=g_tree, inp=f_tree, fmap=gmap)
    assert g_tree.norm() == pytest.approx(g_ref.norm(), rel=epsilon)
    assert g_tree(r_0) == pytest.approx(g_ref(r_0), rel=epsilon)


def test_FunctionMapOperator():
    G = vp.FunctionMap(fmap=gmap, prec=epsilon)
    g_tree = G(f_tree)
    assert g_tree.norm() == pytest.approx(g_ref.norm(), rel=epsilon)
    assert g_tree(r_0) == pytest.approx(g_ref(r_0), rel=epsilon)
