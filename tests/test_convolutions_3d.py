import numpy as np
import vampyr as vp
import pytest

epsilon = 1.0e-3
mu = epsilon/10

D = 3
k = 5
N = -2
legendre = vp.LegendreBasis(order=k)
world = vp.D3.BoundingBox(scale=N)
mra = vp.D3.MultiResolutionAnalysis(box=world, basis=legendre)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
ffunc = vp.D3.GaussFunc(coef=alpha, exp=beta, pos=r0)
ref_energy = ffunc.calcCoulombEnergy(ffunc)

ftree = vp.D3.FunctionTree(mra)
vp.build_grid(out=ftree, inp=ffunc)
vp.project(prec=epsilon, out=ftree, inp=ffunc)

def test_Identity():
    I = vp.D3.IdentityConvolution(mra, prec=epsilon)
    gtree = vp.D3.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

def test_Poisson():
    P = vp.D3.PoissonOperator(mra, prec=epsilon)
    gtree = vp.D3.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=P, inp=ftree)
    assert vp.dot(gtree, ftree) == pytest.approx(ref_energy, rel=epsilon)

def test_Helmholtz():
    H = vp.D3.HelmholtzOperator(mra, exp=mu, prec=epsilon)
    gtree = vp.D3.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=H, inp=ftree)
    assert vp.dot(gtree, ftree) == pytest.approx(ref_energy, rel=epsilon)

def test_PeriodicIdentity():
    world = vp.D3.BoundingBox(pbc=True, corner=[-1,-1,-1], nboxes=[2,2,2])
    pbc = vp.D3.MultiResolutionAnalysis(box=world, basis=legendre)

    pfunc = ffunc.periodify([1.0,1.0,1.0])
    ftree = vp.D3.FunctionTree(mra=pbc)
    vp.build_grid(out=ftree, inp=pfunc)
    vp.project(prec=epsilon, out=ftree, inp=pfunc)

    I1 = vp.D3.IdentityConvolution(mra=pbc, prec=epsilon, root=-5)
    gtree1 = vp.D3.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree1, oper=I1, inp=ftree)
    assert gtree1.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)

    I2 = vp.D3.IdentityConvolution(mra=pbc, prec=epsilon, reach=5)
    gtree2 = vp.D3.FunctionTree(mra=pbc)
    vp.apply(prec=epsilon, out=gtree2, oper=I2, inp=ftree)
    assert gtree2.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)
