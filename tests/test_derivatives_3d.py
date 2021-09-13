import numpy as np
import vampyr as vp
import pytest

epsilon = 1.0e-3

D = 3
k = 5
N = -2
interpol = vp.InterpolatingBasis(order=k)
world = vp.D3.BoundingBox(scale=N)
mra = vp.D3.MultiResolutionAnalysis(box=world, basis=interpol)

r0 = [0.8, 0.8, 0.8]
beta = 10.0
alpha = (beta/np.pi)**(D/2.0)
g = vp.D3.GaussFunc(coef=alpha, exp=beta, pos=r0)
gx = g.differentiate(dir=0)
gy = g.differentiate(dir=1)
gz = g.differentiate(dir=2)

f = vp.D3.FunctionTree(mra)
vp.build_grid(out=f, inp=g)
vp.project(prec=epsilon, out=f, inp=g)

fx = vp.D3.FunctionTree(mra)
vp.build_grid(out=fx, inp=gx)
vp.project(prec=epsilon, out=fx, inp=gx)

fy = vp.D3.FunctionTree(mra)
vp.build_grid(out=fy, inp=gy)
vp.project(prec=epsilon, out=fy, inp=gy)

fz = vp.D3.FunctionTree(mra)
vp.build_grid(out=fz, inp=gz)
vp.project(prec=epsilon, out=fz, inp=gz)

def test_Gradient():
    D = vp.D3.ABGVOperator(mra, a=0.0, b=0.0)
    grad_f = vp.gradient(oper=D, inp=f)

    assert len(grad_f) == 3
    assert grad_f[0].integrate() == pytest.approx(fx.integrate(), abs=epsilon)
    assert grad_f[1].integrate() == pytest.approx(fy.integrate(), abs=epsilon)
    assert grad_f[2].integrate() == pytest.approx(fz.integrate(), abs=epsilon)
    assert grad_f[0].getSquareNorm() == pytest.approx(fx.getSquareNorm(), rel=epsilon)
    assert grad_f[1].getSquareNorm() == pytest.approx(fy.getSquareNorm(), rel=epsilon)
    assert grad_f[2].getSquareNorm() == pytest.approx(fz.getSquareNorm(), rel=epsilon)

def test_Divergence():
    D = vp.D3.ABGVOperator(mra, a=0.0, b=0.0)
    grad_f = []
    grad_f.append(f)
    grad_f.append(f)
    grad_f.append(f)
    lap_f = vp.D3.FunctionTree(mra)
    vp.divergence(out=lap_f, oper=D, inp=grad_f)

    ref_vec = []
    ref_vec.append(fx)
    ref_vec.append(fy)
    ref_vec.append(fz)
    ref_lap = vp.D3.FunctionTree(mra)
    vp.build_grid(out=ref_lap, inp=ref_vec)
    vp.add(out=ref_lap, inp=ref_vec)

    assert lap_f.integrate() == pytest.approx(ref_lap.integrate(), abs=epsilon)
    assert lap_f.getSquareNorm() == pytest.approx(ref_lap.getSquareNorm(), rel=epsilon)
