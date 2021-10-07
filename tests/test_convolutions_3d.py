import numpy as np
from vampyr import vampyr3d as vp
import pytest

epsilon = 1.0e-3
mu = epsilon/10

D = 3
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
ffunc = vp.GaussFunc(coef=alpha, exp=beta, pos=r0)
ref_energy = ffunc.calcCoulombEnergy(ffunc)

ftree = vp.FunctionTree(mra)
vp.build_grid(out=ftree, inp=ffunc)
vp.project(prec=epsilon, out=ftree, inp=ffunc)

def test_Identity():
    I = vp.IdentityConvolution(mra, prec=epsilon)

    gtree = vp.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    gtree2 = I(ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

def test_Poisson():
    P = vp.PoissonOperator(mra, prec=epsilon)

    gtree = vp.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=P, inp=ftree)
    assert vp.dot(gtree, ftree) == pytest.approx(ref_energy, rel=epsilon)

    gtree2 = P(ftree)
    assert vp.dot(gtree2, ftree) == pytest.approx(ref_energy, rel=epsilon)

def test_Helmholtz():
    H = vp.HelmholtzOperator(mra, exp=mu, prec=epsilon)

    gtree = vp.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=H, inp=ftree)
    assert vp.dot(gtree, ftree) == pytest.approx(ref_energy, rel=epsilon)

    gtree2 = H(ftree)
    assert vp.dot(gtree2, ftree) == pytest.approx(ref_energy, rel=epsilon)

def test_PeriodicIdentity():
    world = vp.BoundingBox(pbc=True, corner=[-1,-1,-1], nboxes=[2,2,2])
    pbc = vp.MultiResolutionAnalysis(box=world, order=k)

    pfunc = ffunc.periodify([1.0,1.0,1.0])
    ftree = vp.FunctionTree(mra=pbc)
    vp.build_grid(out=ftree, inp=pfunc)
    vp.project(prec=epsilon, out=ftree, inp=pfunc)

    I1 = vp.IdentityConvolution(mra=pbc, prec=epsilon, root=-5)
    gtree1 = vp.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree1, oper=I1, inp=ftree)
    assert gtree1.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    I2 = vp.IdentityConvolution(mra=pbc, prec=epsilon, reach=5)
    gtree2 = vp.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree2, oper=I2, inp=ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)
