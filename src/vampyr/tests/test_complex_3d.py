import numpy as np
import pytest

from vampyr import LegendreBasis
from vampyr import vampyr3d as vp3

# The tree, arithmetic and projection layers are dimension-generic templates,
# so this file is deliberately a transcription of the dimension-generic half
# of test_complex_1d.py rather than a new set of physics checks. Anything
# genuinely 1D (time evolution, the Schrodinger semigroup) stays there.

r0 = [0.5, 0.5, 0.5]
sigma = 0.01
k0 = [3.0, -2.0, 1.0]

order = 5
precision = 1.0e-5

lbasis = LegendreBasis(order=order)
box = vp3.BoundingBox(0)
mra = vp3.MultiResolutionAnalysis(box, lbasis)

P = vp3.ScalingProjector(mra, prec=precision)
Pc = vp3.ScalingProjector(mra, prec=precision, dtype=complex)


def gauss(r):
    d2 = sum((r[i] - r0[i]) ** 2 for i in range(3))
    return np.exp(-d2 / sigma)


def wavepacket(r):
    phase = sum(k0[i] * r[i] for i in range(3))
    return gauss(r) * np.exp(1j * phase)


psi = Pc(wavepacket)

# Independently projected real/imaginary parts. These are separate adaptive
# projections, so they live on their own grids and agree with psi only to
# within the projection error -- use them where projection agreement is the
# thing under test, never for exactness claims.
u = P(lambda r: wavepacket(r).real)
v = P(lambda r: wavepacket(r).imag)

# Same-grid real/imaginary parts, extracted from psi itself. Use these when
# asserting a binding invariant (complex path == 2x2 real composition), so
# that grid differences cannot masquerade as a binding error.
ur = psi.real()
vi = psi.imag()


def test_dtype_dispatch():
    assert isinstance(Pc, vp3.ComplexScalingProjector)
    assert isinstance(psi, vp3.ComplexFunctionTree)
    with pytest.raises(TypeError):
        vp3.ScalingProjector(mra, prec=precision, dtype=int)


def test_complex_projection_pointwise():
    for r in ([0.5, 0.5, 0.5], [0.45, 0.5, 0.55], [0.55, 0.45, 0.5]):
        assert psi(r) == pytest.approx(wavepacket(r), rel=1e-6, abs=1e-8)


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


def test_from_real_imag_roundtrip():
    """Exactness is only claimed for the same-grid round trip. Assembling from
    independently projected parts agrees with psi to the projection error,
    which at prec=1e-5 is several orders of magnitude looser."""
    exact = vp3.from_real_imag(ur, vi)
    assert (exact - psi).norm() == pytest.approx(0.0, abs=1e-12)
    assert (vp3.ComplexFunctionTree(ur, vi) - psi).norm() == pytest.approx(0.0, abs=1e-12)

    assembled = vp3.from_real_imag(u, v)
    assert (assembled - psi).norm() / psi.norm() == pytest.approx(0.0, abs=1e-5)


def test_mixed_real_complex_arithmetic():
    """Real operands are promoted in place; the result must match doing the
    promotion by hand."""
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
    """conj() is exact on the existing grid, so compare against psi's own
    parts -- not against separately projected ones."""
    c = psi.conj()
    assert (c.real() - ur).norm() == pytest.approx(0.0, abs=1e-12)
    assert (c.imag() + vi).norm() == pytest.approx(0.0, abs=1e-12)
    assert (c.conj() - psi).norm() == pytest.approx(0.0, abs=1e-12)


def test_dot_conventions():
    """dot conjugates the bra for complex trees, so dot(psi, psi) is the
    squared L2 norm and is real and positive."""
    sesq = vp3.dot(psi, psi)
    assert sesq.imag == pytest.approx(0.0, abs=1e-10)
    assert sesq.real == pytest.approx(psi.squaredNorm(), rel=1e-8)
    bilinear = vp3.dot(psi.conj(), psi)
    assert bilinear == pytest.approx((psi * psi).integrate(), rel=1e-6)
    assert vp3.dot(u, psi) == pytest.approx(vp3.dot(u.to_complex(), psi), rel=1e-12)
    assert vp3.dot(psi, u) == pytest.approx(vp3.dot(psi, u.to_complex()), rel=1e-12)


def test_complex_function_map():
    F = vp3.FunctionMap(fmap=lambda z: z * z, prec=precision, dtype=complex)
    assert isinstance(F, vp3.ComplexFunctionMap)
    sq = F(psi)
    ref = psi * psi
    assert (sq - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_function_map_dtype_dispatch():
    """The real path must be untouched by the dtype wrapper."""
    F = vp3.FunctionMap(fmap=lambda x: x * x, prec=precision)
    assert isinstance(F, vp3._RealFunctionMap)
    assert (F(u) - u * u).norm() / (u * u).norm() == pytest.approx(0.0, abs=1e-6)
    with pytest.raises(TypeError):
        vp3.FunctionMap(fmap=lambda x: x, prec=precision, dtype=int)


def test_identity_convolution_complex():
    """A real operator applied to a complex tree must equal the 2x2 real block
    composition. This is the binding invariant, independent of how accurate
    the kernel itself is at this precision."""
    ident = vp3.IdentityConvolution(mra, precision)
    out = ident(psi)
    assert isinstance(out, vp3.ComplexFunctionTree)
    ref = vp3.from_real_imag(ident(ur), ident(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_derivative_complex():
    D = vp3.ABGVDerivative(mra, 0.5, 0.5)
    for direction in range(3):
        dpsi = D(psi, direction)
        ref = vp3.from_real_imag(D(ur, direction), D(vi, direction))
        assert (dpsi - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_poisson_complex():
    """Poisson is real and linear, so applying it to a complex density is the
    same as applying it to the real and imaginary parts separately."""
    Pop = vp3.PoissonOperator(mra, precision)
    out = Pop(psi)
    assert isinstance(out, vp3.ComplexFunctionTree)
    ref = vp3.from_real_imag(Pop(ur), Pop(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)


def test_helmholtz_complex():
    Hop = vp3.HelmholtzOperator(mra, 1.0, precision)
    out = Hop(psi)
    ref = vp3.from_real_imag(Hop(ur), Hop(vi))
    assert (out - ref).norm() / ref.norm() == pytest.approx(0.0, abs=1e-6)
