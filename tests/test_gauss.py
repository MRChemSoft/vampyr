import numpy as np
import vampyr as vp
import pytest

numprec = 1.0e-10

def test_GaussFunc_1D():
    beta = 100.0
    alpha = (np.pi/beta)**(1.0/2.0)
    r_0 = [0.25]
    r_1 = [0.35]
    f = vp.D1.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    ref = alpha*np.exp(-1.0)
    assert f(r_0) == alpha
    assert f(r_1) == pytest.approx(ref, rel=numprec)
    assert f.getCoef() == alpha
    assert f.getPos() == r_0
    assert f.getExp() == beta
    assert f.getPow(dim=0) == 0

def test_GaussFunc_2D():
    beta = 100.0
    alpha = (np.pi/beta)**(2.0/2.0)
    r_0 = [0.25, 0.25]
    r_1 = [0.25, 0.35]
    f = vp.D2.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    ref = alpha*np.exp(-1.0)
    assert f(r_0) == alpha
    assert f(r_1) == pytest.approx(ref, rel=numprec)
    assert f.getCoef() == alpha
    assert f.getPos() == r_0
    assert f.getExp() == beta
    assert f.getPow(dim=1) == 0

def test_GaussFunc_3D():
    beta = 100.0
    alpha = (np.pi/beta)**(3.0/2.0)
    r_0 = [0.25, 0.25, 0.25]
    r_1 = [0.25, 0.25, 0.35]
    f = vp.D3.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    ref = alpha*np.exp(-1.0)
    assert f(r_0) == alpha
    assert f(r_1) == pytest.approx(ref, rel=numprec)
    assert f.getCoef() == alpha
    assert f.getPos() == r_0
    assert f.getExp() == beta
    assert f.getPow(dim=2) == 0

def test_GaussDerivative_1D():
    beta = 100.0
    alpha = (np.pi/beta)**(1.0/2.0)
    r_0 = [0.25]
    r_1 = [0.35]
    f = vp.D1.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    df = f.differentiate(0)
    ref = -(2.0/10.0)*alpha*beta*np.exp(-1.0)
    assert df.getPow(dim=0) == 1
    assert df(r_0) == 0.0
    assert df(r_1) == pytest.approx(ref, rel=numprec)

def test_GaussDerivative_2D():
    beta = 100.0
    alpha = (np.pi/beta)**(2/2)
    r_0 = [0.25, 0.25]
    r_1 = [0.25, 0.35]
    f = vp.D2.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    df = f.differentiate(dim=1)
    ref = -(2/10)*alpha*beta*np.exp(-1.0)
    assert df.getPow(dim=1) == 1
    assert df(r_0) == 0.0
    assert df(r_1) == pytest.approx(ref, rel=numprec)

def test_GaussDerivative_3D():
    beta = 100.0
    alpha = (np.pi/beta)**(3.0/2.0)
    r_0 = [0.25, 0.25, 0.25]
    r_1 = [0.25, 0.25, 0.35]
    f = vp.D3.GaussFunc(coef=alpha, exp=beta, pos=r_0)
    df = f.differentiate(dim=2)
    ref = -(2.0/10.0)*alpha*beta*np.exp(-1.0)
    assert df.getPow(dim=2) == 1
    assert df(r_0) == 0.0
    assert df(r_1) == pytest.approx(ref, rel=numprec)

def test_GaussOverlap_1D():
    b1 = 10.0
    b2 = 20.0
    r1 = [0.00]
    r2 = [0.25]
    f1 = vp.D1.GaussFunc(exp=b1, pos=r1)
    f2 = vp.D1.GaussFunc(exp=b2, pos=r2)
    p = b1 + b2
    mu = b1*b2/p
    R12 = sum([(x1 - x2)**2 for x1, x2 in zip(r1, r2)])
    ref = (np.pi/p)**(1.0/2.0)*np.exp(-mu*R12)
    assert f1.calcOverlap(f2) == pytest.approx(ref, rel=numprec)

def test_GaussOverlap_2D():
    b1 = 10.0
    b2 = 20.0
    r1 = [0.00, 0.00]
    r2 = [0.25, 0.25]
    f1 = vp.D2.GaussFunc(exp=b1, pos=r1)
    f2 = vp.D2.GaussFunc(exp=b2, pos=r2)
    p = b1 + b2
    mu = b1*b2/p
    R12 = sum([(x1 - x2)**2 for x1, x2 in zip(r1, r2)])
    ref = (np.pi/p)**(2.0/2.0)*np.exp(-mu*R12)
    assert f1.calcOverlap(f2) == pytest.approx(ref, rel=numprec)

def test_GaussOverlap_3D():
    b1 = 10.0
    b2 = 20.0
    r1 = [0.00, 0.00, 0.00]
    r2 = [0.25, 0.25, 0.25]
    f1 = vp.D3.GaussFunc(exp=b1, pos=r1)
    f2 = vp.D3.GaussFunc(exp=b2, pos=r2)
    p = b1 + b2
    mu = b1*b2/p
    R12 = sum([(x1 - x2)**2 for x1, x2 in zip(r1, r2)])
    ref = (np.pi/p)**(3.0/2.0)*np.exp(-mu*R12)
    assert f1.calcOverlap(f2) == pytest.approx(ref, rel=numprec)

def test_GaussNorm_1D():
    beta = 10.0
    f = vp.D1.GaussFunc(exp=beta)
    ref = (np.pi/(2.0*beta))**(1.0/2.0)
    assert f.calcSquareNorm() == pytest.approx(ref, rel=numprec)
    assert f.calcSquareNorm() == f.calcOverlap(f)

def test_GaussNorm_2D():
    beta = 10.0
    f = vp.D2.GaussFunc(exp=beta)
    ref = (np.pi/(2.0*beta))**(2.0/2.0)
    assert f.calcSquareNorm() == pytest.approx(ref, rel=numprec)
    assert f.calcSquareNorm() == f.calcOverlap(f)

def test_GaussNorm_3D():
    beta = 10.0
    f = vp.D3.GaussFunc(exp=beta)
    ref = (np.pi/(2.0*beta))**(3.0/2.0)
    assert f.calcSquareNorm() == pytest.approx(ref, rel=numprec)
    assert f.calcSquareNorm() == f.calcOverlap(f)

def test_GaussEnergy():
    beta = 100.0
    alpha = (np.pi/beta)**(3.0/2.0)
    f = vp.D3.GaussFunc(coef=alpha, exp=beta)
    ref = np.sqrt(2.0*beta/np.pi)
    assert f.calcCoulombEnergy(f) == pytest.approx(ref, rel=numprec)

def test_GaussExp_1D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1]
    r1 = [-0.1]
    f0 = vp.D1.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D1.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D1.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0(r0) + f1(r0)
    assert fexp.size() == 2
    assert fexp.getFunc(term=0).getExp() == b0
    assert fexp.getFunc(term=1).getExp() == b1
    assert fexp(r0) == pytest.approx(ref, rel=numprec)

def test_GaussExp_2D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1]
    r1 = [-0.1, -0.1]
    f0 = vp.D2.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D2.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D2.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0(r0) + f1(r0)
    assert fexp.size() == 2
    assert fexp.getFunc(term=0).getExp() == b0
    assert fexp.getFunc(term=1).getExp() == b1
    assert fexp(r0) == pytest.approx(ref, rel=numprec)

def test_GaussExp_3D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    f0 = vp.D3.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D3.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D3.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0(r0) + f1(r0)
    assert fexp.size() == 2
    assert fexp.getFunc(term=0).getExp() == b0
    assert fexp.getFunc(term=1).getExp() == b1
    assert fexp(r0) == pytest.approx(ref, rel=numprec)

def test_GaussExpEnergy():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    f0 = vp.D3.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D3.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D3.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0.calcCoulombEnergy(f0)
    ref += f0.calcCoulombEnergy(f1)
    ref += f1.calcCoulombEnergy(f0)
    ref += f1.calcCoulombEnergy(f1)
    assert fexp.calcCoulombEnergy() == pytest.approx(ref, rel=numprec)

def test_GaussExpNorm_1D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1]
    r1 = [-0.1]
    f0 = vp.D1.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D1.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D1.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0.calcSquareNorm() + 2.0*f0.calcOverlap(f1) + f1.calcSquareNorm()
    assert fexp.calcSquareNorm() == pytest.approx(ref, rel=numprec)

def test_GaussExpNorm_2D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1]
    r1 = [-0.1, -0.1]
    f0 = vp.D2.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D2.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D2.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0.calcSquareNorm() + 2.0*f0.calcOverlap(f1) + f1.calcSquareNorm()
    assert fexp.calcSquareNorm() == pytest.approx(ref, rel=numprec)

def test_GaussExpNorm_3D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    f0 = vp.D3.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D3.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D3.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    ref = f0.calcSquareNorm() + 2.0*f0.calcOverlap(f1) + f1.calcSquareNorm()
    assert fexp.calcSquareNorm() == pytest.approx(ref, rel=numprec)

def test_GaussExpDerivative_1D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1]
    r1 = [-0.1]
    r2 = [0.0]
    f0 = vp.D1.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D1.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D1.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    df0 = f0.differentiate(dim=0)
    df1 = f1.differentiate(dim=0)
    dfexp = fexp.differentiate(dim=0)
    ref = df0(r2) + df1(r2)
    assert dfexp(r2) == pytest.approx(ref, rel=numprec)

def test_GaussExpDerivative_2D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1]
    r1 = [-0.1, -0.1]
    r2 = [0.0, 0.0]
    f0 = vp.D2.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D2.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D2.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    df0 = f0.differentiate(dim=1)
    df1 = f1.differentiate(dim=1)
    dfexp = fexp.differentiate(dim=1)
    ref = df0(r2) + df1(r2)
    assert dfexp(r2) == pytest.approx(ref, rel=numprec)

def test_GaussExpDerivative_3D():
    b0 = 10.0
    b1 = 20.0
    r0 = [0.1, 0.1, 0.1]
    r1 = [-0.1, -0.1, -0.1]
    r2 = [0.0, 0.0, 0.0]
    f0 = vp.D3.GaussFunc(exp=b0, pos=r0)
    f1 = vp.D3.GaussFunc(exp=b1, pos=r1)
    fexp = vp.D3.GaussExp()
    fexp.append(f0)
    fexp.append(f1)
    df0 = f0.differentiate(dim=2)
    df1 = f1.differentiate(dim=2)
    dfexp = fexp.differentiate(dim=2)
    ref = df0(r2) + df1(r2)
    assert dfexp(r2) == pytest.approx(ref, rel=numprec)
