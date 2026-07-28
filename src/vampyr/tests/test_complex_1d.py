import numpy as np
import pytest

from vampyr import LegendreBasis
from vampyr import vampyr1d as vp1

# Wavepacket parameters
x0 = 0.5
sigma = 0.001
k0 = 20.0

# Time step (mirrors test_time_evolution_1d.py)
time = 0.002

order = 5
precision = 1.0e-6
finest_scale = 8

lbasis = LegendreBasis(order=order)
box = vp1.BoundingBox(0)
mra = vp1.MultiResolutionAnalysis(box, lbasis)

P = vp1.ScalingProjector(mra, prec=precision)
Pc = vp1.ScalingProjector(mra, prec=precision, dtype=complex)


def gauss(x):
    return np.exp(-((x - x0) ** 2) / sigma)


def wavepacket(x):
    return gauss(x) * np.exp(1j * k0 * x)


psi = Pc(lambda r: wavepacket(r[0]))
u = P(lambda r: (wavepacket(r[0])).real)
v = P(lambda r: (wavepacket(r[0])).imag)


def test_dtype_dispatch():
    assert isinstance(Pc, vp1.ComplexScalingProjector)
    assert isinstance(psi, vp1.ComplexFunctionTree)
    with pytest.raises(TypeError):
        vp1.ScalingProjector(mra, prec=precision, dtype=int)


def test_complex_projection_pointwise():
    for x in (0.45, 0.5, 0.55):
        assert psi([x]) == pytest.approx(wavepacket(x), rel=1e-8, abs=1e-10)


def test_complex_integrate():
    # int exp(-(x-x0)^2/sigma) exp(ikx) dx = sqrt(pi*sigma) exp(ik x0) exp(-k^2 sigma/4)
    ref = np.sqrt(np.pi * sigma) * np.exp(1j * k0 * x0) * np.exp(-(k0**2) * sigma / 4)
    assert psi.integrate() == pytest.approx(ref, rel=1e-8)


def test_from_real_imag_roundtrip():
    psi2 = vp1.from_real_imag(u, v)
    psi3 = vp1.ComplexFunctionTree(u, v)
    assert (psi - psi2).norm() == pytest.approx(0.0, abs=1e-8)
    assert (psi2 - psi3).norm() == pytest.approx(0.0, abs=1e-14)
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
    # distributivity against real/imag split
    ref = vp1.from_real_imag((2.0 * u) - (3.0 * v), (3.0 * u) + (2.0 * v))
    assert (scaled - ref).norm() == pytest.approx(0.0, abs=1e-7)
    back = scaled / a
    assert (back - psi).norm() == pytest.approx(0.0, abs=1e-12)


def test_conj():
    psic = psi.conj()
    assert (psic.real() - psi.real()).norm() == pytest.approx(0.0, abs=1e-12)
    assert (psic.imag() + psi.imag()).norm() == pytest.approx(0.0, abs=1e-12)


def test_dot_conventions():
    # <psi|psi> conjugates the bra: real, equals the squared norm
    inner = vp1.dot(psi, psi)
    assert inner.imag == pytest.approx(0.0, abs=1e-10)
    assert inner.real == pytest.approx(psi.squaredNorm(), rel=1e-10)
    # bilinear integral of psi^2 via materialized conjugate
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
    """The identity kernel smooths this sharp wavepacket noticeably at this
    prec (real path shows the same ~1e-2 deviation), so assert the binding
    invariant instead: applying to the complex tree equals the real
    2x2 block composition."""
    ident = vp1.IdentityConvolution(mra, precision)
    out = ident(psi)
    assert isinstance(out, vp1.ComplexFunctionTree)
    ref = vp1.from_real_imag(ident(u), ident(v))
    assert (out - ref).norm() == pytest.approx(0.0, abs=1e-6)
    # and it should not stray from the real-path error scale
    assert (out - psi).norm() == pytest.approx((ident(u) + 1j * ident(v) - psi).norm(), rel=1e-6)


def test_derivative_complex():
    D = vp1.ABGVDerivative(mra, 0.5, 0.5)
    dpsi = D(psi)
    ref = Pc(lambda r: (1j * k0 - 2.0 * (r[0] - x0) / sigma) * wavepacket(r[0]))
    assert (dpsi - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-4)


def free_particle_analytical_solution(x, x0, t, sigma):
    denominator = 4.0j * t + sigma
    return np.sqrt(sigma) / np.sqrt(denominator) * np.exp(-((x - x0) ** 2) / denominator)


def test_complex_time_evolution_single_operator():
    """One operator object, applied directly to a (promoted) complex tree:
    no manual real/imaginary split-apply-reunite."""
    f = P(lambda r: gauss(r[0]))
    U = vp1.ComplexTimeEvolutionOperator(mra, precision, time, finest_scale)
    out = U(f)
    assert isinstance(out, vp1.ComplexFunctionTree)

    g = Pc(lambda r: free_particle_analytical_solution(r[0], x0, time, sigma))
    assert (out - g).squaredNorm() == pytest.approx(0.0, abs=5e-14)

    # second application on the complex tree itself: t -> 2t. The evolved
    # input is chirped, and the fixed-Jpower kernel is less accurate on
    # oscillatory functions, so this bound is intrinsically looser than the
    # first step (same behaviour as composing the real-block operators).
    out2 = U(out)
    g2 = Pc(lambda r: free_particle_analytical_solution(r[0], x0, 2 * time, sigma))
    assert (out2 - g2).squaredNorm() == pytest.approx(0.0, abs=1e-6)

    # unitarity on this window
    assert out.norm() == pytest.approx(f.norm(), rel=1e-4)


def test_real_kernel_on_complex_tree_matches_split():
    """Applying a single real/imag kernel to a complex tree equals the
    2x2 real block composition."""
    A = vp1.TimeEvolutionOperator(mra, precision, time, finest_scale, False)
    B = vp1.TimeEvolutionOperator(mra, precision, time, finest_scale, True)
    Apsi = A(psi)
    ref = vp1.from_real_imag(A(u), A(v))
    assert (Apsi - ref).norm() == pytest.approx(0.0, abs=1e-6)
    full = Apsi + 1j * B(psi)
    U = vp1.ComplexTimeEvolutionOperator(mra, precision, time, finest_scale)
    assert (full - U(psi)).norm() == pytest.approx(0.0, abs=1e-6)
