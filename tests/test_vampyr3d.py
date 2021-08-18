import numpy as np
import vampyr as vp
import pytest


order = 5
prec = 1e-3
min_scale = -4
max_depth = 25
corner = [-1, -1, -1]
boxes = [2, 2, 2]

print("MRCPP VERSION : " + vp.mrcpp_version())

# basis = vp.InterpolatingBasis(order)
# world = vp.D3.BoundingBox(min_scale, corner, boxes)
# MRA = vp.D3.MultiResolutionAnalysis(world, basis, max_depth)

#def d_phi_exact(x):
#    beta = 100.0
#    alpha = (beta/np.pi)**(3/2)
#
#    return -2.0*beta*alpha*x[0]*np.exp(-beta*(x[0]**2 + x[1]**2 + x[2]**2))

# H = vp.D3.HelmholtzOperator(MRA, mu, prec)
# P = vp.D3.PoissonOperator(MRA, prec)

# inp_tree_p = vp.D3.FunctionTree(MRA)
# out_tree_p = vp.D3.FunctionTree(MRA)

# def f_exact(r):
#     R2 = (r[0]-r_0[0])**2 + (r[1]-r_0[1])**2 + (r[2]-r_0[2])**2
#     return alpha*np.exp(-beta*R2)
#
# beta = 1000.0
# alpha = (beta/np.pi)**(3.0/2.0)
# r_0 = [0.0, 0.0, np.pi/2.0]
# f_func = vp.D3.GaussFunc(beta, alpha, r_0)
#
# f_tree = vp.D3.FunctionTree(MRA)
# vp.build_grid(f_tree, f_func)
# vp.project(prec, f_tree, f_exact)
#
# P = vp.D3.PoissonOperator(MRA, prec)
# g_tree = vp.D3.FunctionTree(MRA)
# vp.apply(prec, g_tree, P, f_tree)

# tree_generator(inp_tree_p, MRA, prec, beta, mid, mu)
# vp.apply(prec, out_tree_p, P, inp_tree_p)

# print(inp_tree_p)
# print(out_tree_p)
#
# integral = inp_tree_p.integrate()
# print("integral: {}".format(integral))
#
# energy = vp.dot(out_tree_p, inp_tree_p)
# print("energy: {}".format(energy))
#
# val = out_tree_p.evalf([0, 0, 0])
# ref = poisson_exact([0, 0, 0])
# print("f(0,0,0):  {}".format(val))
# print("reference: {}".format(ref))

#v_generator(v_tree_pois, MRA, prec, beta, mid)

#add_tree = vp.vampyr3d.FunctionTree(MRA)
#add_vec_tree = vp.vampyr3d.FunctionTree(MRA)
#mult_vec_tree = vp.vampyr3d.FunctionTree(MRA)

#vp.apply(prec, phi_tree, H, v_tree)


# def test_BBGetScale():
#     assert world.getScale() == min_scale
#
# def test_IBGetScalingOrder():
#     assert basis.getScalingOrder() == order
#
# def test_MRAGetOrder():
#     assert MRA.getOrder() == order
#
# def test_GaussEvalf():
#     r_1 = [np.pi, np.pi, np.pi]
#     assert f_func(r_0) == pytest.approx(f_exact(r_0))
#     assert f_func(r_1) == pytest.approx(f_exact(r_1))
#
# def test_TreeEvalf():
#     r_1 = [np.pi, np.pi, np.pi]
#     assert f_tree.evalf(r_0) == pytest.approx(f_exact(r_0), rel=prec)
#     assert f_tree.evalf(r_1) == pytest.approx(f_exact(r_1), rel=prec)
#
# def test_TreeIntegral():
#     assert f_tree.integrate() == pytest.approx(1.0, rel=prec)
#
# def test_TreeDotProduct():
#     ref_energy = f_func.calcCoulombEnergy(f_func)
#     num_energy = vp.dot(f_tree, g_tree)
#     assert num_energy == pytest.approx(ref_energy, rel=prec)

# def test_evalf_helm():
#     assert phi_tree.evalf([0, 0, 0]) == pytest.approx(phi_exact([0, 0, 0]), rel=prec)
#
# def test_add():
#     vp.add(prec/10, add_tree, 1.0, phi_tree, -1, phi_tree_pois)
#     assert add_tree.evalf([0, 0, 0]) == pytest.approx(0.0, abs=prec*10)
#
#
# def test_add_vec():
#     sum_vec = []
#     sum_vec.append(tuple([1.0, phi_tree]))
#     sum_vec.append(tuple([-1.0, phi_tree]))
#     vp.add(prec/10, add_vec_tree, sum_vec)
#     assert add_vec_tree.evalf([0.0, 0.0, 0.0]) == pytest.approx(0.0, abs=prec*10)
#
#
# def test_multiply():
#     mult_tree = vp.vampyr3d.FunctionTree(MRA)
#     vp.multiply(prec, mult_tree, 1, phi_tree, phi_tree_pois)
#     assert mult_tree.evalf([0, 0, 0]) == pytest.approx(
#                    phi_exact([0, 0, 0])**2, rel=prec)
#
#
# def test_multiply_vec():
#     multiply_vec = []
#     multiply_vec.append(tuple([1.0, phi_tree]))
#     multiply_vec.append(tuple([1.0, phi_tree]))
#     vp.multiply(prec/10, mult_vec_tree, multiply_vec)
#     assert mult_vec_tree.evalf([0, 0, 0]) == pytest.approx(
#                    phi_exact([0, 0, 0])**2, rel=prec)
#
#
# def test_gaussFunc():
#     sigma = 0.01
#     pos = [0.0, 0.0, 0.0]
#     power = [0, 0, 0]
#     alpha = 1.0/(2.0*np.pi*sigma**2)**(3/2)
#     beta = 1.0/(2.0*sigma**2)
#
#     gauss = vp.vampyr3d.GaussFunc(beta, alpha, pos, power)
#     g_tree = vp.vampyr3d.FunctionTree(MRA)
#     vp.build_grid(g_tree, gauss)
#     vp.project(prec, g_tree, gauss)
#     assert g_tree.integrate() == pytest.approx(1.0, rel=prec)
#
#
# def test_gaussFuncAnalytic():
#     def gauss(x):
#         beta = 100.0
#         return np.exp(-beta * (x[0]**2 + x[1]**2 + x[2]**2))
#
#     g_tree = vp.vampyr3d.FunctionTree(MRA)
#     vp.project(prec, g_tree, gauss)
#     assert g_tree.integrate() == pytest.approx(0.00556832799694068, rel=prec)
#
#
# def test_gaussExp():
#     gf = vp.vampyr3d.GaussFunc(1.0, 1.0)
#     gexp = vp.vampyr3d.GaussExp()
#     gexp.append(gf)
#
#     assert gexp.evalf([0.0, 0.0, 0.0]) == pytest.approx(1.0, rel=prec)
#
#
# def test_power_square_and_dot():
#     tmp_1_tree = vp.vampyr3d.FunctionTree(MRA)
#     tmp_2_tree = vp.vampyr3d.FunctionTree(MRA)
#     tmp_3_tree = vp.vampyr3d.FunctionTree(MRA)
#
#     vp.power(prec, tmp_1_tree, phi_tree, 2.0)
#     vp.square(prec, tmp_2_tree, phi_tree)
#     assert vp.dot(phi_tree, phi_tree) == pytest.approx(tmp_2_tree.integrate())
#     assert tmp_1_tree.integrate() == pytest.approx(tmp_2_tree.integrate(),
#                    rel=prec)
#
#
# def test_vec_dot():
#     vec = []
#     vec.append(tuple([1.0, phi_tree]))
#     vec.append(tuple([-1.0, phi_tree]))
#     tmp_1_tree = vp.vampyr3d.FunctionTree(MRA)
#     tmp_2_tree = vp.vampyr3d.FunctionTree(MRA)
#     vp.square(prec, tmp_2_tree, phi_tree)
#     vp.dot(prec, tmp_1_tree, vec, vec)
#     assert tmp_1_tree.integrate() == pytest.approx(tmp_2_tree.integrate()*2.0,
#                    rel=prec)
#
#
# def test_gradient_and_divergence():
#     grad_tree = vp.vampyr3d.FunctionTree(MRA)
#     out_grad_tree = vp.vampyr3d.FunctionTree(MRA)
#     pos = [0.0, 0.0, 0.0]
#     power = [0, 0, 0]
#     beta = 100.0
#     alpha = (beta/np.pi)**(3/2)
#     gauss = vp.vampyr3d.GaussFunc(beta, alpha, pos, power)
#     vp.project(prec, grad_tree, gauss)
#     D = vp.vampyr3d.ABGVOperator(MRA, 0.0, 0.0)
#     grad = vp.gradient(D, grad_tree)
#     assert grad[0][1].evalf([0.1, 0.0, 0.0]) == pytest.approx(
#                    d_phi_exact([0.1, 0.0, 0.0]), rel=prec)
#     grad_tree_vec = []
#     grad_tree_vec.append(tuple([1.0, grad_tree]))
#     grad_tree_vec.append(tuple([1.0, grad_tree]))
#     grad_tree_vec.append(tuple([1.0, grad_tree]))
#     vp.divergence(out_grad_tree, D, grad_tree_vec)
#     assert out_grad_tree.evalf([0.1, 0.1, 0.1]) == pytest.approx(
#                    3.0*grad[0][1].evalf([0.1, 0.1, 0.1]), rel=prec)
#
#
# def test_copy_func():
#     copy_tree = vp.vampyr3d.FunctionTree(MRA)
#     vp.copy_grid(copy_tree, phi_tree)
#     vp.copy_func(copy_tree, phi_tree)
#     assert copy_tree.integrate() == pytest.approx(phi_tree.integrate(), rel=prec)


## FIXME This test isn't really testing anything
#def test_function_tree_squared():
#    power = [2, 2, 2]
#    f = vp.GaussFunc(beta, 4.0*beta**2, [mid, mid, mid], power)
#
#    f_tree = vp.FunctionTree(MRA)
#    vp.project(prec, f_tree, f)
#    f_tree.evalf([0.0, 0.0, 0.0])
#    f_tree.square()
#    assert f_tree.evalf([0.0, 0.0, 0.0]) == pytest.approx(0.0)
#
#
## FIXME This test isn't really testing anything
#def test_function_tree_power():
#    power = [0, 0, 0]
#    f = vp.GaussFunc(beta, 4.0, [mid, mid, mid], power)
#
#    f_tree = vp.FunctionTree(MRA)
#    vp.project(prec, f_tree, f)
#    f_tree.evalf([0.0, 0.0, 0.0])
#    f_tree.power(2.0)
#    assert f_tree.evalf([0.0, 0.0, 0.0]) == pytest.approx(16.0)
#
#
## FIXME This test isn't really testing anything
#def test_function_tree_add():
#    power = [0, 0, 0]
#    f = vp.GaussFunc(beta, 4.0, [mid, mid, mid], power)
#
#    f_tree = vp.FunctionTree(MRA)
#    vp.project(prec, f_tree, f)
#    f_tree.evalf([0.0, 0.0, 0.0])
#    f_tree.add(2.0, f_tree)
#    assert f_tree.evalf([0.0, 0.0, 0.0]) == pytest.approx(12.0)
#
#
## FIXME This test isn't really testing anything
#def test_function_tree_multiply():
#    power = [0, 0, 0]
#    f = vp.GaussFunc(beta, 4.0, [mid, mid, mid], power)
#
#    f_tree = vp.FunctionTree(MRA)
#    vp.project(prec, f_tree, f)
#    f_tree.evalf([0.0, 0.0, 0.0])
#    f_tree.multiply(1.0, f_tree)
#    assert f_tree.evalf([0.0, 0.0, 0.0]) == pytest.approx(16.0)
