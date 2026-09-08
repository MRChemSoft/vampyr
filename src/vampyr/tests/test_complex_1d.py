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


def test_time_evolution_single_operator():
    """One complex kernel, one convolution, no real/imaginary split."""
    f = P(lambda r: gauss(r[0]))
    U = vp1.TimeEvolutionOperator(mra, precision, time, finest_scale)
    assert U.iscomplex()

    out = U(f)
    assert isinstance(out, vp1.ComplexFunctionTree)

    g = Pc(lambda r: free_particle_analytical_solution(r[0], x0, time, sigma))
    assert (out - g).squaredNorm() == pytest.approx(0.0, abs=5e-14)

    # t -> 2t. The evolved input is chirped and the fixed-Jpower kernel is
    # less accurate on oscillatory functions, hence the looser bound.
    out2 = U(out)
    g2 = Pc(lambda r: free_particle_analytical_solution(r[0], x0, 2 * time, sigma))
    assert (out2 - g2).squaredNorm() == pytest.approx(0.0, abs=1e-6)

    # unitarity on this window
    assert out.norm() == pytest.approx(f.norm(), rel=1e-4)


def test_real_input_is_promoted():
    """MRCPP aborts on a real tree, so the binding promotes first. The result
    must match promoting by hand."""
    f_real = P(lambda r: gauss(r[0]))
    f_cplx = Pc(lambda r: complex(gauss(r[0]), 0.0))

    U = vp1.TimeEvolutionOperator(mra, precision, time, finest_scale)
    from_real = U(f_real)
    from_cplx = U(f_cplx)

    assert isinstance(from_real, vp1.ComplexFunctionTree)
    assert (from_real - from_cplx).norm() / from_cplx.norm() == pytest.approx(0.0, abs=1e-10)

    # a real state picks up an imaginary part
    assert from_real.imag().squaredNorm() > 1e-12


def test_adaptive_and_uniform_agree():
    """The adaptive build is the default now, and must land on the same
    semigroup as the fixed-scale one."""
    adaptive = vp1.TimeEvolutionOperator(mra, precision, time)
    uniform = vp1.TimeEvolutionOperator(mra, precision, time, finest_scale)
    assert adaptive.iscomplex()
    assert vp1.TimeEvolutionOperator.Adaptive == -1

    a = adaptive(psi)
    b = uniform(psi)
    assert (a - b).norm() / b.norm() == pytest.approx(0.0, abs=1e-6)


def test_max_Jpower_reaches_the_adaptive_path():
    """max_Jpower is honoured with or without a finest_scale, so a truncated
    expansion must move the answer."""
    full = vp1.TimeEvolutionOperator(mra, precision, time, max_Jpower=30)
    short = vp1.TimeEvolutionOperator(mra, precision, time, max_Jpower=2)
    assert (full(psi) - short(psi)).norm() > 1e-10


def test_legacy_imaginary_flag_is_rejected():
    """Both old signatures must raise rather than binding the flag to
    max_Jpower, which the integer caster would otherwise accept."""
    with pytest.raises(TypeError):
        vp1.TimeEvolutionOperator(mra, precision, time, True)
    with pytest.raises(TypeError):
        vp1.TimeEvolutionOperator(mra, precision, time, finest_scale, False)


def test_real_convolutions_still_real():
    """Only the time evolution kernel went complex. The rest keep a real
    expansion, and a real tree in still gives a real tree out."""
    I = vp1.IdentityConvolution(mra, prec=precision)
    assert I.isreal()
    f = P(lambda r: gauss(r[0]))
    assert isinstance(I(f), vp1.FunctionTree)
    assert isinstance(I(psi), vp1.ComplexFunctionTree)


def test_grid_helpers():
    """build_grid, refine_grid and clear_grid on complex trees."""
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
    """square(conjugate=True) is the density, imaginary part is round-off."""
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
    """1D divergence of a one-component field is just the derivative."""
    D = vp1.ABGVDerivative(mra, 0.0, 0.0)
    div = vp1.advanced.divergence(oper=D, inp=[psi])

    d_psi = vp1.ComplexFunctionTree(mra)
    vp1.advanced.apply(out=d_psi, oper=D, inp=psi, dir=0)

    assert (div - d_psi).norm() / d_psi.norm() == pytest.approx(0.0, abs=1e-10)


def test_time_evolution_rejects_unsupported_mra():
    """A negative root scale would index past the power integrals and
    dereference a null pointer in MRCPP."""
    shifted = vp1.MultiResolutionAnalysis(vp1.BoundingBox(-2), lbasis)
    with pytest.raises(ValueError):
        vp1.TimeEvolutionOperator(shifted, precision, time)
