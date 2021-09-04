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
beta = 100.0
alpha = (beta/np.pi)**(D/2.0)
ffunc = vp.D1.GaussFunc(coef=alpha, exp=beta, pos=r0)

ftree = vp.D1.FunctionTree(mra)
vp.build_grid(out=ftree, inp=ffunc)
vp.project(prec=epsilon, out=ftree, inp=ffunc)

def test_Identity():
    I = vp.D1.IdentityConvolution(mra, prec=epsilon)
    gtree = vp.D1.FunctionTree(mra)
    vp.apply(prec=epsilon, out=gtree, oper=I, inp=ftree)
    assert gtree.integrate() == pytest.approx(ftree.integrate(), rel=epsilon)
