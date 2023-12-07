import numpy as np
import pytest

from vampyr import vampyr1d as vp1
from vampyr import vampyr3d as vp

epsilon = 1.0e-3
mu = epsilon / 10

D = 3
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta / np.pi) ** (D / 2.0)
ffunc = vp.GaussFunc(alpha=alpha, beta=beta, position=r0)
ref_energy = ffunc.calcCoulombEnergy(ffunc)

ftree = vp.FunctionTree(mra)
vp.advanced.build_grid(out=ftree, inp=ffunc)
vp.advanced.project(prec=epsilon, out=ftree, inp=ffunc)


def test_GaussKernel():
    b = 1.0e4
    a = (b / np.pi) ** (D / 2.0)
    ifunc = vp1.GaussFunc(alpha=a, beta=b)
    iexp = vp1.GaussExp()
    iexp.append(ifunc)
    I = vp.ConvolutionOperator(mra, iexp, prec=epsilon)

    ftree = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=ftree, inp=ffunc)
    vp.advanced.project(prec=epsilon, out=ftree, inp=ffunc)

    gtree = vp.FunctionTree(mra)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    gtree2 = I(ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)


def test_CartesianConvolution():
    b = 1.0e4
    a = (b / np.pi) ** (D / 2.0)
    ifunc = vp1.GaussFunc(alpha=a, beta=b)
    iexp = vp1.GaussExp()
    iexp.append(ifunc)
    O = vp.CartesianConvolution(mra, iexp, prec=epsilon)
    O.setCartesianComponents(0, 0, 0)

    ftree = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=ftree, inp=ffunc)
    vp.advanced.project(prec=epsilon, out=ftree, inp=ffunc)

    gtree = vp.FunctionTree(mra)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=O, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    gtree2 = O(ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)


def test_Identity():
    I = vp.IdentityConvolution(mra, prec=epsilon)

    gtree = vp.FunctionTree(mra)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    gtree2 = I(ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)


def test_Poisson():
    P = vp.PoissonOperator(mra, prec=epsilon)

    gtree = vp.FunctionTree(mra)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=P, inp=ftree)
    assert vp.dot(gtree, ftree) == pytest.approx(ref_energy, rel=epsilon)

    gtree2 = P(4 * np.pi * ftree)
    assert vp.dot(gtree2, ftree) == pytest.approx(ref_energy, rel=epsilon)


def test_Helmholtz():
    H = vp.HelmholtzOperator(mra, exp=mu, prec=epsilon)

    gtree = vp.FunctionTree(mra)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=H, inp=ftree)
    assert vp.dot(gtree, ftree) == pytest.approx(ref_energy, rel=epsilon)

    gtree2 = H(4 * np.pi * ftree)
    assert vp.dot(gtree2, ftree) == pytest.approx(ref_energy, rel=epsilon)


def test_PeriodicIdentity():
    world = vp.BoundingBox(pbc=True, corner=[-1, -1, -1], nboxes=[2, 2, 2])
    pbc = vp.MultiResolutionAnalysis(box=world, order=k)

    pfunc = ffunc.periodify([1.0, 1.0, 1.0])
    ftree = vp.FunctionTree(mra=pbc)
    vp.advanced.build_grid(out=ftree, inp=pfunc)
    vp.advanced.project(prec=epsilon, out=ftree, inp=pfunc)

    I1 = vp.IdentityConvolution(mra=pbc, prec=epsilon, root=-5)
    gtree1 = vp.FunctionTree(mra=pbc)
    vp.advanced.apply(prec=epsilon, out=gtree1, oper=I1, inp=ftree)
    assert gtree1.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    I2 = vp.IdentityConvolution(mra=pbc, prec=epsilon, reach=5)
    gtree2 = vp.FunctionTree(mra=pbc)
    vp.advanced.apply(prec=epsilon, out=gtree2, oper=I2, inp=ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)
