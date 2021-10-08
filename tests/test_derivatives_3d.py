import numpy as np
import pytest

from vampyr import vampyr3d as vp

epsilon = 1.0e-3

D = 3
k = 5
N = -2
world = vp.BoundingBox(scale=N)
mra = vp.MultiResolutionAnalysis(box=world, order=k)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta / np.pi) ** (D / 2.0)
g = vp.GaussFunc(coef=alpha, exp=beta, pos=r0)
gx = g.differentiate(dir=0)
gy = g.differentiate(dir=1)
gz = g.differentiate(dir=2)

f = vp.FunctionTree(mra)
vp.build_grid(out=f, inp=g)
vp.project(prec=epsilon, out=f, inp=g)

fx = vp.FunctionTree(mra)
vp.build_grid(out=fx, inp=gx)
vp.project(prec=epsilon, out=fx, inp=gx)

fy = vp.FunctionTree(mra)
vp.build_grid(out=fy, inp=gy)
vp.project(prec=epsilon, out=fy, inp=gy)

fz = vp.FunctionTree(mra)
vp.build_grid(out=fz, inp=gz)
vp.project(prec=epsilon, out=fz, inp=gz)


def test_Gradient():
    D = vp.ABGVDerivative(mra, a=0.0, b=0.0)
    grad_f = vp.gradient(oper=D, inp=f)

    assert len(grad_f) == 3
    assert grad_f[0].integrate() == pytest.approx(fx.integrate(), abs=epsilon)
    assert grad_f[1].integrate() == pytest.approx(fy.integrate(), abs=epsilon)
    assert grad_f[2].integrate() == pytest.approx(fz.integrate(), abs=epsilon)
    assert grad_f[0].norm() == pytest.approx(fx.norm(), rel=epsilon)
    assert grad_f[1].norm() == pytest.approx(fy.norm(), rel=epsilon)
    assert grad_f[2].norm() == pytest.approx(fz.norm(), rel=epsilon)


def test_Divergence():
    D = vp.ABGVDerivative(mra, a=0.0, b=0.0)
    grad_f = []
    grad_f.append(f)
    grad_f.append(f)
    grad_f.append(f)
    lap_f = vp.divergence(oper=D, inp=grad_f)

    ref_vec = []
    ref_vec.append(fx)
    ref_vec.append(fy)
    ref_vec.append(fz)
    ref_lap = vp.FunctionTree(mra)
    vp.build_grid(out=ref_lap, inp=ref_vec)
    vp.add(out=ref_lap, inp=ref_vec)

    assert lap_f.integrate() == pytest.approx(ref_lap.integrate(), abs=epsilon)
    assert lap_f.norm() == pytest.approx(ref_lap.norm(), rel=epsilon)


def test_OverloadedOperators():
    D = vp.ABGVDerivative(mra, a=0.0, b=0.0)
    assert D.getOrder() == 1

    gx = D(f, axis=0)
    gy = D(f, axis=1)
    gz = D(f, axis=2)
    assert gx.integrate() == pytest.approx(fx.integrate(), abs=epsilon)
    assert gy.integrate() == pytest.approx(fy.integrate(), abs=epsilon)
    assert gz.integrate() == pytest.approx(fz.integrate(), abs=epsilon)
    assert gx.norm() == pytest.approx(fx.norm(), rel=epsilon)
    assert gy.norm() == pytest.approx(fy.norm(), rel=epsilon)
    assert gz.norm() == pytest.approx(fz.norm(), rel=epsilon)
