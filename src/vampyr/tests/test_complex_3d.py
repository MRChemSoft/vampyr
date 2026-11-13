import numpy as np
import pytest

from vampyr import LegendreBasis
from vampyr import vampyr1d as vp1
from vampyr import vampyr3d as vp3

# The dimension-generic part of test_complex_1d.py. Time evolution in 3D is in
# test_time_evolution_3d.py.

r0 = [0.5, 0.5, 0.5]
sigma = 0.01
k0 = [3.0, -2.0, 1.0]

order = 5
precision = 1.0e-5

lbasis = LegendreBasis(order=order)
box = vp3.BoundingBox(0)
mra = vp3.MultiResolutionAnalysis(box, lbasis)

P = vp3.ScalingProjector(mra, prec=precision)
Pc = vp3.ComplexScalingProjector(mra, prec=precision)


def gauss(r):
    d2 = sum((r[i] - r0[i]) ** 2 for i in range(3))
    return np.exp(-d2 / sigma)


def wavepacket(r):
    phase = sum(k0[i] * r[i] for i in range(3))
    return gauss(r) * np.exp(1j * phase)


psi = Pc(wavepacket)

# u and v are projected on their own grids, so they equal the parts of psi only
# to the projection error.
u = P(lambda r: wavepacket(r).real)
v = P(lambda r: wavepacket(r).imag)

# The parts of psi itself, on its grid.
ur = psi.real()
vi = psi.imag()


def test_complex_projector():
    assert isinstance(psi, vp3.ComplexFunctionTree)
    # A complex tree is projected directly, not point by point through Python
    again = Pc(psi)
    assert (again - psi).norm() / psi.norm() == pytest.approx(0.0, abs=1e-5)


def test_complex_projection_pointwise():
    for r in ([0.5, 0.5, 0.5], [0.45, 0.5, 0.55], [0.55, 0.45, 0.5]):
        assert psi(r) == pytest.approx(wavepacket(r), rel=1e-6, abs=1e-8)


def test_real_and_imag_norms():
    assert ur.squaredNorm() + vi.squaredNorm() == pytest.approx(
        psi.squaredNorm(), rel=1e-12
    )


def test_complex_integrate():
    # separable: the 3D integral is the product of three 1D Gaussian integrals
    ref = 1.0 + 0.0j
    for i in range(3):
        ref *= (
            np.sqrt(np.pi * sigma)
            * np.exp(1j * k0[i] * r0[i])
            * np.exp(-(k0[i] ** 2) * sigma / 4)
        )
    assert psi.integrate() == pytest.approx(ref, rel=1e-6)


def test_from_real_and_imaginary_parts():
    """Exact for the parts of psi; for separately projected parts, accurate to
    the projection error."""
    exact = vp3.ComplexFunctionTree(ur, vi)
    assert (exact - psi).norm() == pytest.approx(0.0, abs=1e-12)

    assembled = vp3.ComplexFunctionTree(u, v)
    assert (assembled - psi).norm() / psi.norm() == pytest.approx(0.0, abs=1e-5)


def test_mixed_real_complex_arithmetic():
    """A real operand is promoted; the result matches promoting it by hand."""
    uc = u.to_complex()
    assert ((psi + u) - (psi + uc)).norm() == pytest.approx(0.0, abs=1e-12)
    assert ((psi - u) - (psi - uc)).norm() == pytest.approx(0.0, abs=1e-12)
    assert ((u + psi) - (uc + psi)).norm() == pytest.approx(0.0, abs=1e-12)
    assert ((psi * u) - (psi * uc)).norm() == pytest.approx(0.0, abs=1e-12)
    assert ((u * psi) - (uc * psi)).norm() == pytest.approx(0.0, abs=1e-12)


def test_complex_scalar_ops():
    assert ((1j * psi) - (psi * 1j)).norm() == pytest.approx(0.0, abs=1e-12)
    assert ((psi * 2.0) - (psi + psi)).norm() == pytest.approx(0.0, abs=1e-10)
    assert ((psi / 2.0) * 2.0 - psi).norm() == pytest.approx(0.0, abs=1e-10)
    # multiplying by i twice is a sign flip
    assert ((1j * (1j * psi)) + psi).norm() == pytest.approx(0.0, abs=1e-10)


def test_conj():
    """conj() keeps the grid, so compare with the parts of psi itself."""
    c = psi.conj()
    assert (c.real() - ur).norm() == pytest.approx(0.0, abs=1e-12)
    assert (c.imag() + vi).norm() == pytest.approx(0.0, abs=1e-12)
    assert (c.conj() - psi).norm() == pytest.approx(0.0, abs=1e-12)


def test_dot_conventions():
    """dot conjugates the bra, so dot(psi, psi) is the squared norm."""
    sesq = vp3.dot(psi, psi)
    assert sesq.imag == pytest.approx(0.0, abs=1e-10)
    assert sesq.real == pytest.approx(psi.squaredNorm(), rel=1e-8)
    bilinear = vp3.dot(psi.conj(), psi)
    assert bilinear == pytest.approx((psi * psi).integrate(), rel=1e-6)
    assert vp3.dot(u, psi) == pytest.approx(vp3.dot(u.to_complex(), psi), rel=1e-12)
    assert vp3.dot(psi, u) == pytest.approx(vp3.dot(psi, u.to_complex()), rel=1e-12)


def test_complex_function_map():
    F = vp3.ComplexFunctionMap(fmap=lambda z: z * z, prec=precision)
    sq = F(psi)
    ref = psi * psi
    assert (sq - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_identity_convolution_complex():
    """A real operator acts on the real and imaginary parts separately."""
    ident = vp3.IdentityConvolution(mra, precision)
    out = ident(psi)
    assert isinstance(out, vp3.ComplexFunctionTree)
    ref = vp3.ComplexFunctionTree(ident(ur), ident(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_derivative_complex():
    D = vp3.ABGVDerivative(mra, 0.5, 0.5)
    for direction in range(3):
        dpsi = D(psi, direction)
        ref = vp3.ComplexFunctionTree(D(ur, direction), D(vi, direction))
        assert (dpsi - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_poisson_complex():
    """Poisson is real, so it acts on the real and imaginary parts separately."""
    Pop = vp3.PoissonOperator(mra, precision)
    out = Pop(psi)
    assert isinstance(out, vp3.ComplexFunctionTree)
    ref = vp3.ComplexFunctionTree(Pop(ur), Pop(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_helmholtz_complex():
    Hop = vp3.HelmholtzOperator(mra, 1.0, precision)
    out = Hop(psi)
    ref = vp3.ComplexFunctionTree(Hop(ur), Hop(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_cartesian_convolution_complex():
    """CartesianConvolution has its own __call__, so it needs its own complex overload."""
    beta = 1.0e4
    alpha = (beta / np.pi) ** 1.5
    kernel = vp1.GaussExp()
    kernel.append(vp1.GaussFunc(alpha=alpha, beta=beta))
    O = vp3.CartesianConvolution(mra, kernel, prec=precision)
    O.setCartesianComponents(0, 0, 0)
    out = O(psi)
    assert isinstance(out, vp3.ComplexFunctionTree)
    ref = vp3.ComplexFunctionTree(O(ur), O(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)
