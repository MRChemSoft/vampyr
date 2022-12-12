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
ffunc = vp.GaussFunc(alpha=alpha, beta=beta, position=r0)


def test_Identity():
    I = vp.IdentityConvolution(mra, prec=epsilon)

    ftree = vp.FunctionTree(mra)
    vp.advanced.build_grid(out=ftree, inp=ffunc)
    vp.advanced.project(prec=epsilon, out=ftree, inp=ffunc)

    gtree = vp.FunctionTree(mra)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    gtree2 = I(ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)


def test_PeriodicIdentity():
    world = vp.BoundingBox(pbc=True)
    pbc = vp.MultiResolutionAnalysis(box=world, order=k)

    pfunc = ffunc.periodify([1.0])
    ftree = vp.FunctionTree(mra=pbc)
    vp.advanced.build_grid(out=ftree, inp=pfunc)
    vp.advanced.project(prec=epsilon, out=ftree, inp=pfunc)

    I = vp.IdentityConvolution(mra=pbc, prec=epsilon, root=0, reach=5)
    gtree = vp.FunctionTree(mra=pbc)
    vp.advanced.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)


def test_PeriodicIdentity():
    world = vp.BoundingBox(pbc=True, corner=[-1], nboxes=[2])
    pbc = vp.MultiResolutionAnalysis(box=world, order=k)

    pfunc = ffunc.periodify([1.0])
    ftree = vp.FunctionTree(mra=pbc)
    vp.advanced.build_grid(out=ftree, inp=pfunc)
    vp.advanced.project(prec=epsilon, out=ftree, inp=pfunc)

    I1 = vp.IdentityConvolution(mra=pbc, prec=epsilon, root=-2)
    gtree1 = vp.FunctionTree(mra=pbc)
    vp.advanced.apply(prec=epsilon, out=gtree1, oper=I1, inp=ftree)
    assert gtree1.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    I2 = vp.IdentityConvolution(mra=pbc, prec=epsilon, reach=5)
    gtree2 = vp.FunctionTree(mra=pbc)
    vp.advanced.apply(prec=epsilon, out=gtree2, oper=I2, inp=ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)
