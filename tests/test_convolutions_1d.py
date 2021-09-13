import numpy as np
import vampyr as vp
import pytest

epsilon = 1.0e-3

D = 1
k = 5
N = -2
legendre = vp.InterpolatingBasis(order=k)
world = vp.D1.BoundingBox(scale=N)
mra = vp.D1.MultiResolutionAnalysis(box=world, basis=legendre)

r0 = [0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
ffunc = vp.D1.GaussFunc(coef=alpha, exp=beta, pos=r0)

def test_Identity():
    ftree = vp.D1.FunctionTree(mra)
    vp.build_grid(out=ftree, inp=ffunc)
    vp.project(prec=epsilon, out=ftree, inp=ffunc)

    I = vp.D1.IdentityConvolution(mra, prec=epsilon)
    gtree = vp.D1.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

def test_PeriodicIdentity():
    world = vp.D1.BoundingBox(pbc=True)
    pbc = vp.D1.MultiResolutionAnalysis(box=world, basis=legendre)

    pfunc = ffunc.periodify([1.0])
    ftree = vp.D1.FunctionTree(mra=pbc)
    vp.build_grid(out=ftree, inp=pfunc)
    vp.project(prec=epsilon, out=ftree, inp=pfunc)

    I = vp.D1.IdentityConvolution(mra=pbc, prec=epsilon, root=0, reach=5)
    gtree = vp.D1.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

def test_PeriodicIdentity():
    world = vp.D1.BoundingBox(pbc=True, corner=[-1], nboxes=[2])
    pbc = vp.D1.MultiResolutionAnalysis(box=world, basis=legendre)

    pfunc = ffunc.periodify([1.0])
    ftree = vp.D1.FunctionTree(mra=pbc)
    vp.build_grid(out=ftree, inp=pfunc)
    vp.project(prec=epsilon, out=ftree, inp=pfunc)

    I1 = vp.D1.IdentityConvolution(mra=pbc, prec=epsilon, root=-2)
    gtree1 = vp.D1.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree1, oper=I1, inp=ftree)
    assert gtree1.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    I2 = vp.D1.IdentityConvolution(mra=pbc, prec=epsilon, reach=5)
    gtree2 = vp.D1.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree2, oper=I2, inp=ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)
