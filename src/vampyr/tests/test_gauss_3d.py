import numpy as np
import pytest

from vampyr import vampyr3d as vp

numprec = 1.0e-10


def test_GaussFunc():
    beta = 100.0
    alpha = (np.pi / beta) ** (3.0 / 2.0)
    r_0 = [0.25, 0.25, 0.25]
    r_1 = [0.25, 0.25, 0.35]
    f = vp.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    ref = alpha * np.exp(-1.0)
    assert f(r_0) == alpha
    assert f(r_1) == pytest.approx(ref, rel=numprec)
    assert f.coef() == alpha
    assert f.pos() == r_0
    assert f.exp() == beta
    assert f.pow(dim=2) == 0


def test_GaussDerivative():
    beta = 100.0
    alpha = (np.pi / beta) ** (3.0 / 2.0)
    r_0 = [0.25, 0.25, 0.25]
    r_1 = [0.25, 0.25, 0.35]
    f = vp.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    df = f.differentiate(dir=2)
    ref = -(2.0 / 10.0) * alpha * beta * np.exp(-1.0)
    assert df(r_0) == 0.0
    assert df(r_1) == pytest.approx(ref, rel=numprec)


def test_GaussOverlap():
    b1 = 10.0
    b2 = 20.0
    r1 = [0.00, 0.00, 0.00]
    r2 = [0.25, 0.25, 0.25]
    f1 = vp.GaussFunc(exp=b1, pos=r1)
    f2 = vp.GaussFunc(exp=b2, pos=r2)
    p = b1 + b2
    mu = b1 * b2 / p
    R12 = sum([(x1 - x2) ** 2 for x1, x2 in zip(r1, r2)])
    ref = (np.pi / p) ** (3.0 / 2.0) * np.exp(-mu * R12)
    assert f1.overlap(f2) == pytest.approx(ref, rel=numprec)


def test_GaussNorm():
    beta = 10.0
    f = vp.GaussFunc(exp=beta)
    ref = (np.pi / (2.0 * beta)) ** (3.0 / 2.0)
    assert f.squaredNorm() == pytest.approx(ref, rel=numprec)
    assert f.squaredNorm() == f.overlap(f)


def test_GaussEnergy():
    beta = 100.0
    alpha = (np.pi / beta) ** (3.0 / 2.0)
    f = vp.GaussFunc(coef=alpha, exp=beta)
    ref = np.sqrt(2.0 * beta / np.pi)
    assert f.calcCoulombEnergy(f) == pytest.approx(ref, rel=numprec)


def test_GaussExp():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    f0 = vp.GaussFunc(exp=b0, pos=r0)
    f1 = vp.GaussFunc(exp=b1, pos=r1)
    fexp = vp.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0(r0) + f1(r0)
    assert fexp.size() == 2
    assert fexp.func(term=0).exp() == b0
    assert fexp.func(term=1).exp() == b1
    assert fexp(r0) == pytest.approx(ref, rel=numprec)


def test_GaussExpEnergy():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    f0 = vp.GaussFunc(exp=b0, pos=r0)
    f1 = vp.GaussFunc(exp=b1, pos=r1)
    fexp = vp.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0.calcCoulombEnergy(f0)
    ref += f0.calcCoulombEnergy(f1)
    ref += f1.calcCoulombEnergy(f0)
    ref += f1.calcCoulombEnergy(f1)
    assert fexp.calcCoulombEnergy() == pytest.approx(ref, rel=numprec)


def test_GaussExpNorm():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    f0 = vp.GaussFunc(exp=b0, pos=r0)
    f1 = vp.GaussFunc(exp=b1, pos=r1)
    fexp = vp.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0.squaredNorm() + 2.0 * f0.overlap(f1) + f1.squaredNorm()
    assert fexp.squaredNorm() == pytest.approx(ref, rel=numprec)


def test_GaussExpDerivative():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    r2 = [0.0, 0.0, 0.0]
    f0 = vp.GaussFunc(exp=b0, pos=r0)
    f1 = vp.GaussFunc(exp=b1, pos=r1)
    fexp = vp.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    df0 = f0.differentiate(dir=2)
    df1 = f1.differentiate(dir=2)
    dfexp = fexp.differentiate(dir=2)
    ref = df0(r2) + df1(r2)
    assert dfexp(r2) == pytest.approx(ref, rel=numprec)
