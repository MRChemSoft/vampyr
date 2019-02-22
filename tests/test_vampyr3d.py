import os
import sys
sys.path.append(os.getenv('VAMPYR_MODULE_PATH'))

from math import isclose

from pot_gen import v_generator
import numpy as np
import vampyr3d as vp


min_scale = -4
max_depth = 25
order = 5
prec = 1e-3
beta = 100.0
mu = 10.0
corner = [-1, -1, -1]
boxes = [2, 2, 2]

world = vp.BoundingBox(min_scale, corner, boxes)
basis = vp.InterpolatingBasis(order)
MRA = vp.MultiResolutionAnalysis(world, basis, max_depth)

mid = 0.0


def phi_exact(x):
    beta = 100.0
    alpha = (beta/np.pi)**(3/2)

    return alpha*np.exp(-beta*(x[0]**2 + x[1]**2 + x[2]**2))


def d_phi_exact(x):
    beta = 100.0
    alpha = (beta/np.pi)**(3/2)

    return -2.0*beta*alpha*x[0]*np.exp(-beta*(x[0]**2 + x[1]**2 + x[2]**2))


H = vp.HelmholtzOperator(MRA, mu, prec)
P = vp.PoissonOperator(MRA, prec)

phi_tree = vp.FunctionTree(MRA)
phi_tree_pois = vp.FunctionTree(MRA)
v_tree = vp.FunctionTree(MRA)
v_generator(v_tree, MRA, prec, beta, mid, mu)
v_tree_pois = vp.FunctionTree(MRA)
v_generator(v_tree_pois, MRA, prec, beta, mid)

add_tree = vp.FunctionTree(MRA)
add_vec_tree = vp.FunctionTree(MRA)
mult_vec_tree = vp.FunctionTree(MRA)

vp.apply(prec, phi_tree, H, v_tree)
vp.apply(prec, phi_tree_pois, P, v_tree_pois)


def test_IsIntWorking():
    assert isclose(1.0, phi_tree.integrate(), rel_tol=prec)


def test_BBGetScale():
    assert world.getScale() == min_scale


def test_IBGetScalingOrder():
    assert basis.getScalingOrder() == order


def test_MRAGetOrder():
    assert MRA.getOrder() == order


def test_evalf_helm():
    assert isclose(phi_tree.evalf([0, 0, 0]), phi_exact([0, 0, 0]),
                   rel_tol=prec)


def test_evalf_pelm():
    assert isclose(phi_tree_pois.evalf([0, 0, 0]),
                   phi_exact([0, 0, 0]), rel_tol=prec)


def test_add():
    vp.add(prec/10, add_tree, 1.0, phi_tree, -1, phi_tree_pois)
    assert isclose(add_tree.evalf([0, 0, 0]), 0.0, abs_tol=prec*10)


def test_add_vec():
    sum_vec = []
    sum_vec.append(tuple([1.0, phi_tree]))
    sum_vec.append(tuple([-1.0, phi_tree]))
    vp.add(prec/10, add_vec_tree, sum_vec)
    assert isclose(add_vec_tree.evalf([0.0, 0.0, 0.0]), 0.0, abs_tol=prec*10)


def test_multiply():
    mult_tree = vp.FunctionTree(MRA)
    vp.multiply(prec, mult_tree, 1, phi_tree, phi_tree_pois)
    assert isclose(mult_tree.evalf([0, 0, 0]),
                   phi_exact([0, 0, 0])**2, rel_tol=prec)


def test_multiply_vec():
    multiply_vec = []
    multiply_vec.append(tuple([1.0, phi_tree]))
    multiply_vec.append(tuple([1.0, phi_tree]))
    vp.multiply(prec/10, mult_vec_tree, multiply_vec)
    assert isclose(mult_vec_tree.evalf([0, 0, 0]),
                   phi_exact([0, 0, 0])**2, rel_tol=prec)


def test_gaussFunc():
    sigma = 0.01
    pos = [0.0, 0.0, 0.0]
    power = [0, 0, 0]
    alpha = 1.0/(2.0*np.pi*sigma**2)**(3/2)
    beta = 1.0/(2.0*sigma**2)

    gauss = vp.GaussFunc(beta, alpha, pos, power)
    g_tree = vp.FunctionTree(MRA)
    vp.build_grid(g_tree, gauss)
    vp.project(prec, g_tree, gauss)
    assert isclose(g_tree.integrate(), 1.0, rel_tol=prec)


def test_gaussExp():
    gf = vp.GaussFunc(1.0, 1.0)
    gexp = vp.GaussExp()
    gexp.append(gf)

    assert isclose(gexp.evalf([0.0, 0.0, 0.0]), 1.0, rel_tol=prec)


def test_power_square_and_dot():
    tmp_1_tree = vp.FunctionTree(MRA)
    tmp_2_tree = vp.FunctionTree(MRA)
    tmp_3_tree = vp.FunctionTree(MRA)

    vp.power(prec, tmp_1_tree, phi_tree, 2.0)
    vp.square(prec, tmp_2_tree, phi_tree)
    assert isclose(vp.dot(phi_tree, phi_tree), tmp_2_tree.integrate())
    assert isclose(tmp_1_tree.integrate(), tmp_2_tree.integrate(),
                   rel_tol=prec)


def test_vec_dot():
    vec = []
    vec.append(tuple([1.0, phi_tree]))
    vec.append(tuple([-1.0, phi_tree]))
    tmp_1_tree = vp.FunctionTree(MRA)
    tmp_2_tree = vp.FunctionTree(MRA)
    vp.square(prec, tmp_2_tree, phi_tree)
    vp.dot(prec, tmp_1_tree, vec, vec)
    assert isclose(tmp_1_tree.integrate(), tmp_2_tree.integrate()*2.0,
                   rel_tol=prec)


def test_gradient_and_divergence():
    grad_tree = vp.FunctionTree(MRA)
    out_grad_tree = vp.FunctionTree(MRA)
    vp.project(prec, grad_tree, phi_exact)
    D = vp.ABGVOperator(MRA, 0.0, 0.0)
    grad = vp.gradient(D, grad_tree)
    assert isclose(grad[0][1].evalf([0.1, 0.0, 0.0]),
                   d_phi_exact([0.1, 0.0, 0.0]), rel_tol=prec)
    grad_tree_vec = []
    grad_tree_vec.append(tuple([1.0, grad_tree]))
    grad_tree_vec.append(tuple([1.0, grad_tree]))
    grad_tree_vec.append(tuple([1.0, grad_tree]))
    vp.divergence(out_grad_tree, D, grad_tree_vec)
    assert isclose(out_grad_tree.evalf([0.1, 0.1, 0.1]),
                   3.0*grad[0][1].evalf([0.1, 0.1, 0.1]), rel_tol=prec)


def test_copy_func():
    copy_tree = vp.FunctionTree(MRA)
    vp.copy_grid(copy_tree, phi_tree)
    vp.copy_func(copy_tree, phi_tree)
    assert isclose(copy_tree.integrate(), phi_tree.integrate(), rel_tol=prec)
