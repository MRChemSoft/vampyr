import numpy as np
import vampyr1d as vp

from math import isclose

min_scale = -4
max_depth = 25
order = 7
prec = 1e-5

corner = -1
boxes = 2

world = vp.BoundingBox(min_scale, corner, boxes)

basis = vp.InterpolatingBasis(order)

MRA = vp.MultiResolutionAnalysis(world, basis, max_depth)


def phi(x):
    beta = 100
    alpha = (beta/np.pi)**(1/2)

    return alpha*np.exp(-beta*x**2)

def d_phi(x):
    beta = 100
    alpha = (beta/np.pi)**(1/2)

    return -2.0*alpha*beta*x*np.exp(-beta*x**2)



phi_tree = vp.FunctionTree(MRA)
d_phi_tree = vp.FunctionTree(MRA)

add_tree = vp.FunctionTree(MRA)
mult_tree = vp.FunctionTree(MRA)
diff_tree = vp.FunctionTree(MRA)

vp.project(prec, phi_tree, phi)
vp.project(prec, d_phi_tree, d_phi)

def test_IsIntWorking():
    assert isclose(1.0, phi_tree.integrate(), rel_tol=prec)


def test_BBGetScale():
    assert world.getScale() == min_scale


def test_IBGetScalingOrder():
    assert basis.getScalingOrder() == order


def test_MRAGetOrder():
    assert MRA.getOrder() == order


def test_add():
    vp.add(prec/10, add_tree, 1.0, phi_tree, -1, phi_tree)
    assert isclose(add_tree.evalf(0), 0.0, abs_tol=prec*10)


def test_multiply():
    vp.multiply(prec, mult_tree, 1, phi_tree, phi_tree)
    assert isclose(mult_tree.evalf(0), phi(0)**2, rel_tol=prec)

def test_differentiation():
    D = vp.ABGVOperator(MRA, 0.0, 0.0)
    vp.apply(diff_tree, D, phi_tree, 0)
    assert isclose(diff_tree.evalf(0.2), d_phi_tree.evalf(0.2), rel_tol=prec)
