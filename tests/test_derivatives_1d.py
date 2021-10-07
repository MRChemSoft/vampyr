import numpy as np
from vampyr import vampyr1d as vp
import pytest

epsilon = 1.0e-3

D = 1
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.4]
beta = 100.0
alpha = (beta/np.pi)**(D/2.0)
g = vp.GaussFunc(coef=alpha, exp=beta, pos=r0)
dg = g.differentiate(dir=0)
ddg = dg.differentiate(dir=0)

f = vp.FunctionTree(mra)
vp.build_grid(out=f, inp=g)
vp.project(prec=epsilon/10, out=f, inp=g)

df = vp.FunctionTree(mra)
vp.build_grid(out=df, inp=dg)
vp.project(prec=epsilon/10, out=df, inp=dg)

ddf = vp.FunctionTree(mra)
vp.build_grid(out=ddf, inp=ddg)
vp.project(prec=epsilon/10, out=ddf, inp=ddg)

def test_DerivativeABGV_00():
    D = vp.ABGVDerivative(mra, a=0.0, b=0.0)
    assert D.getOrder() == 1

    dh = vp.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.norm() == pytest.approx(df.norm(), rel=epsilon)

    dh2 = D(f)
    assert dh2.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh2.norm() == pytest.approx(df.norm(), rel=epsilon)

def test_DerivativeABGV_55():
    D = vp.ABGVDerivative(mra, a=0.5, b=0.5)
    assert D.getOrder() == 1

    dh = vp.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.norm() == pytest.approx(df.norm(), rel=epsilon)

    dh2 = D(f)
    assert dh2.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh2.norm() == pytest.approx(df.norm(), rel=epsilon)

def test_DerivativePH_1():
    D = vp.PHDerivative(mra, order=1)
    assert D.getOrder() == 1

    dh = vp.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.norm() == pytest.approx(df.norm(), rel=epsilon)

    dh2 = D(f)
    assert dh2.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh2.norm() == pytest.approx(df.norm(), rel=epsilon)

def test_DerivativePH_2():
    D = vp.PHDerivative(mra, order=2)
    assert D.getOrder() == 2

    ddh = vp.FunctionTree(mra)
    vp.apply(out=ddh, oper=D, inp=f, dir=0)
    assert ddh.integrate() == pytest.approx(ddf.integrate(), abs=epsilon)
    assert ddh.norm() == pytest.approx(ddf.norm(), rel=epsilon)

    ddh2 = D(f)
    assert ddh2.integrate() == pytest.approx(ddf.integrate(), abs=epsilon)
    assert ddh2.norm() == pytest.approx(ddf.norm(), rel=epsilon)

def test_DerivativeBS_1():
    D = vp.BSDerivative(mra, order=1)
    assert D.getOrder() == 1

    dh = vp.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.norm() == pytest.approx(df.norm(), rel=epsilon)

    dh2 = D(f)
    assert dh2.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh2.norm() == pytest.approx(df.norm(), rel=epsilon)

def test_DerivativeBS_2():
    D = vp.BSDerivative(mra, order=2)
    assert D.getOrder() == 2

    ddh = vp.FunctionTree(mra)
    vp.apply(out=ddh, oper=D, inp=f, dir=0)
    assert ddh.integrate() == pytest.approx(ddf.integrate(), abs=epsilon)
    assert ddh.norm() == pytest.approx(ddf.norm(), rel=epsilon)

    ddh2 = D(f)
    assert ddh2.integrate() == pytest.approx(ddf.integrate(), abs=epsilon)
    assert ddh2.norm() == pytest.approx(ddf.norm(), rel=epsilon)
