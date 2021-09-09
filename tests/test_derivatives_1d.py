import numpy as np
import vampyr as vp
import pytest

epsilon = 1.0e-3

D = 1
k = 5
N = -2
interpol = vp.InterpolatingBasis(order=k)
world = vp.D1.BoundingBox(scale=N)
mra = vp.D1.MultiResolutionAnalysis(box=world, basis=interpol)

r0 = [0.4]
beta = 100.0
alpha = (beta/np.pi)**(D/2.0)
g = vp.D1.GaussFunc(coef=alpha, exp=beta, pos=r0)
dg = g.differentiate(dir=0)
ddg = dg.differentiate(dir=0)

f = vp.D1.FunctionTree(mra)
vp.build_grid(out=f, inp=g)
vp.project(prec=epsilon/10, out=f, inp=g)

df = vp.D1.FunctionTree(mra)
vp.build_grid(out=df, inp=dg)
vp.project(prec=epsilon/10, out=df, inp=dg)

ddf = vp.D1.FunctionTree(mra)
vp.build_grid(out=ddf, inp=ddg)
vp.project(prec=epsilon/10, out=ddf, inp=ddg)

def test_DerivativeABGV_00():
    D = vp.D1.ABGVOperator(mra, a=0.0, b=0.0)
    assert D.getOrder() == 1

    dh = vp.D1.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.getSquareNorm() == pytest.approx(df.getSquareNorm(), rel=epsilon)

def test_DerivativeABGV_55():
    D = vp.D1.ABGVOperator(mra, a=0.5, b=0.5)
    assert D.getOrder() == 1

    dh = vp.D1.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.getSquareNorm() == pytest.approx(df.getSquareNorm(), rel=epsilon)

def test_DerivativePH_1():
    D = vp.D1.PHOperator(mra, order=1)
    assert D.getOrder() == 1

    dh = vp.D1.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.getSquareNorm() == pytest.approx(df.getSquareNorm(), rel=epsilon)

def test_DerivativePH_2():
    D = vp.D1.PHOperator(mra, order=2)
    assert D.getOrder() == 2

    ddh = vp.D1.FunctionTree(mra)
    vp.apply(out=ddh, oper=D, inp=f, dir=0)
    assert ddh.integrate() == pytest.approx(ddf.integrate(), abs=epsilon)
    assert ddh.getSquareNorm() == pytest.approx(ddf.getSquareNorm(), rel=epsilon)

def test_DerivativeBS_1():
    D = vp.D1.BSOperator(mra, order=1)
    assert D.getOrder() == 1

    dh = vp.D1.FunctionTree(mra)
    vp.apply(out=dh, oper=D, inp=f, dir=0)
    assert dh.integrate() == pytest.approx(df.integrate(), abs=epsilon)
    assert dh.getSquareNorm() == pytest.approx(df.getSquareNorm(), rel=epsilon)

def test_DerivativeBS_2():
    D = vp.D1.BSOperator(mra, order=2)
    assert D.getOrder() == 2

    ddh = vp.D1.FunctionTree(mra)
    vp.apply(out=ddh, oper=D, inp=f, dir=0)
    assert ddh.integrate() == pytest.approx(ddf.integrate(), abs=epsilon)
    assert ddh.getSquareNorm() == pytest.approx(ddf.getSquareNorm(), rel=epsilon)
