import numpy as np
import pytest

from vampyr import LegendreBasis
from vampyr import vampyr1d as vp1

# Wavepacket parameters
x0 = 0.5
sigma = 0.001
k0 = 20.0

order = 5
precision = 1.0e-6

lbasis = LegendreBasis(order=order)
box = vp1.BoundingBox(0)
mra = vp1.MultiResolutionAnalysis(box, lbasis)

P = vp1.ScalingProjector(mra, prec=precision)
Pc = vp1.ComplexScalingProjector(mra, prec=precision)


def gauss(x):
    return np.exp(-((x - x0) ** 2) / sigma)


def wavepacket(x):
    return gauss(x) * np.exp(1j * k0 * x)


psi = Pc(lambda r: wavepacket(r[0]))
u = P(lambda r: wavepacket(r[0]).real)
v = P(lambda r: wavepacket(r[0]).imag)


def test_complex_projector():
    assert isinstance(psi, vp1.ComplexFunctionTree)
    # A complex tree is projected directly, as the real projector does with a
    # real tree, not point by point through Python.
    again = Pc(psi)
    assert isinstance(again, vp1.ComplexFunctionTree)
    assert (again - psi).norm() == pytest.approx(0.0, abs=1e-6)
    out = vp1.ComplexFunctionTree(mra)
    vp1.advanced.project(prec=precision, out=out, inp=psi)
    assert (out - psi).norm() == pytest.approx(0.0, abs=1e-6)


def test_complex_projection_pointwise():
    for x in (0.45, 0.5, 0.55):
        assert psi([x]) == pytest.approx(wavepacket(x), rel=1e-8, abs=1e-10)


def test_complex_integrate():
    # int exp(-(x-x0)^2/sigma) exp(ikx) dx = sqrt(pi*sigma) exp(ik x0) exp(-k^2 sigma/4)
    ref = np.sqrt(np.pi * sigma) * np.exp(1j * k0 * x0) * np.exp(-(k0**2) * sigma / 4)
    assert psi.integrate() == pytest.approx(ref, rel=1e-8)


def test_from_real_and_imaginary_parts():
    psi2 = vp1.ComplexFunctionTree(u, v)
    assert (psi - psi2).norm() == pytest.approx(0.0, abs=1e-8)
    assert (psi.real() - u).norm() == pytest.approx(0.0, abs=1e-8)
    assert (psi.imag() - v).norm() == pytest.approx(0.0, abs=1e-8)


def test_mixed_real_complex_arithmetic():
    # 1j * (real tree) promotes; real + complex resolves through __radd__
    psi2 = u + 1j * v
    assert isinstance(psi2, vp1.ComplexFunctionTree)
    assert (psi - psi2).norm() == pytest.approx(0.0, abs=1e-8)

    psi3 = u.to_complex() + 1j * v
    assert (psi2 - psi3).norm() == pytest.approx(0.0, abs=1e-13)

    # subtraction from the real side
    diff = u - psi2
    ref = -1j * v
    assert (diff - ref).norm() == pytest.approx(0.0, abs=1e-13)

    # real tree times complex tree
    prod = u * psi
    ref2 = u.to_complex() * psi
    assert (prod - ref2).norm() == pytest.approx(0.0, abs=1e-12)


def test_complex_scalar_ops():
    a = 2.0 + 3.0j
    scaled = a * psi
    assert scaled.norm() == pytest.approx(abs(a) * psi.norm(), rel=1e-12)
    # (2 + 3i)(u + iv) = (2u - 3v) + i(3u + 2v)
    ref = vp1.ComplexFunctionTree((2.0 * u) - (3.0 * v), (3.0 * u) + (2.0 * v))
    assert (scaled - ref).norm() == pytest.approx(0.0, abs=1e-7)
    back = scaled / a
    assert (back - psi).norm() == pytest.approx(0.0, abs=1e-12)


def test_conj():
    psic = psi.conj()
    assert (psic.real() - psi.real()).norm() == pytest.approx(0.0, abs=1e-12)
    assert (psic.imag() + psi.imag()).norm() == pytest.approx(0.0, abs=1e-12)


def test_real_and_imag_norms():
    real, imag = psi.real(), psi.imag()
    assert real.squaredNorm() + imag.squaredNorm() == pytest.approx(
        psi.squaredNorm(), rel=1e-12
    )
    imag.normalize()
    assert imag.squaredNorm() == pytest.approx(1.0, rel=1e-12)


def test_dot_conventions():
    # dot conjugates the bra, so dot(psi, psi) is real and equals the squared norm
    inner = vp1.dot(psi, psi)
    assert inner.imag == pytest.approx(0.0, abs=1e-10)
    assert inner.real == pytest.approx(psi.squaredNorm(), rel=1e-10)
    # without the conjugation: the integral of psi^2
    bilinear = vp1.dot(psi.conj(), psi)
    ref = (psi * psi).integrate()
    assert bilinear == pytest.approx(ref, rel=1e-8)
    # mixed real/complex overloads agree with promotion
    assert vp1.dot(u, psi) == pytest.approx(vp1.dot(u.to_complex(), psi), rel=1e-12)
    assert vp1.dot(psi, u) == pytest.approx(vp1.dot(psi, u.to_complex()), rel=1e-12)


def test_complex_function_map():
    F = vp1.ComplexFunctionMap(fmap=lambda z: z * z, prec=precision)
    sq = F(psi)
    ref = psi * psi
    assert (sq - ref).norm() == pytest.approx(0.0, abs=1e-7)


def test_identity_convolution_complex():
    """A real operator acts on the real and imaginary parts separately."""
    ident = vp1.IdentityConvolution(mra, precision)
    out = ident(psi)
    assert isinstance(out, vp1.ComplexFunctionTree)
    ref = vp1.ComplexFunctionTree(ident(u), ident(v))
    assert (out - ref).norm() == pytest.approx(0.0, abs=1e-6)


def test_derivative_complex():
    D = vp1.ABGVDerivative(mra, 0.5, 0.5)
    dpsi = D(psi)
    ref = Pc(lambda r: (1j * k0 - 2.0 * (r[0] - x0) / sigma) * wavepacket(r[0]))
    assert (dpsi - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-4)


def test_real_convolutions_still_real():
    """The other convolutions are real: a real tree in gives a real tree out."""
    I = vp1.IdentityConvolution(mra, prec=precision)
    assert I.isreal()
    f = P(lambda r: gauss(r[0]))
    assert isinstance(I(f), vp1.FunctionTree)
    assert isinstance(I(psi), vp1.ComplexFunctionTree)


def test_grid_helpers():
    """build_grid, refine_grid and clear_grid on complex trees."""
    out = vp1.ComplexFunctionTree(mra)
    vp1.build_grid(out=out, inp=psi)
    assert out.nNodes() == psi.nNodes()

    out = vp1.ComplexFunctionTree(mra)
    vp1.advanced.build_grid(out=out, inp=psi)
    assert out.nNodes() == psi.nNodes()

    vp1.advanced.copy_func(out=out, inp=psi)
    assert (out - psi).norm() == pytest.approx(0.0, abs=1e-14)

    before = out.nNodes()
    vp1.advanced.refine_grid(out=out, scales=1)
    assert out.nNodes() > before

    # clear_grid drops the coefficients but keeps the grid
    after = out.nNodes()
    vp1.advanced.clear_grid(out=out)
    assert out.nNodes() == after


def test_square_and_power():
    """square(conjugate=True) is the density; its imaginary part is round-off."""
    rho = vp1.ComplexFunctionTree(mra)
    vp1.advanced.square(prec=precision, out=rho, inp=psi, conjugate=True)

    assert rho.imag().norm() == pytest.approx(0.0, abs=1e-10)
    # the integral of |psi|^2 is the squared L2 norm
    assert rho.real().integrate() == pytest.approx(psi.squaredNorm(), rel=1e-6)

    # without the flag it is the algebraic square
    sq = vp1.ComplexFunctionTree(mra)
    vp1.advanced.square(prec=precision, out=sq, inp=psi)
    assert (sq - psi * psi).norm() / sq.norm() == pytest.approx(0.0, abs=1e-6)

    # and __pow__ agrees with square on exponent 2
    assert ((psi**2.0) - sq).norm() / sq.norm() == pytest.approx(0.0, abs=1e-6)


def test_divergence():
    """In 1D the divergence of a one-component field is its derivative."""
    D = vp1.ABGVDerivative(mra, 0.0, 0.0)
    div = vp1.divergence(oper=D, inp=[psi])

    d_psi = vp1.ComplexFunctionTree(mra)
    vp1.advanced.apply(out=d_psi, oper=D, inp=psi, dir=0)

    assert (div - d_psi).norm() / d_psi.norm() == pytest.approx(0.0, abs=1e-10)


def test_sum_and_prod():
    assert (vp1.sum([psi, u.to_complex()]) - (psi + u)).norm() == pytest.approx(
        0.0, abs=1e-12
    )
    assert (
        vp1.sum([(2.0, psi), (1j, psi)]) - (2.0 + 1j) * psi
    ).norm() == pytest.approx(0.0, abs=1e-12)
    assert (vp1.prod([psi, psi]) - psi * psi).norm() == pytest.approx(0.0, abs=1e-12)
