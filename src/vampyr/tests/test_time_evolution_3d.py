import cmath

import pytest

from vampyr import LegendreBasis
from vampyr import vampyr3d as vp3

# One operator tree acts in every direction, so a product of Gaussians evolves
# into the product of the 1D solutions.

order = 5
lbasis = LegendreBasis(order=order)


def free_particle_analytical_solution(x, x0, t, sigma):
    denominator = 4.0j * t + sigma
    return cmath.sqrt(sigma / denominator) * cmath.exp(-((x - x0) ** 2) / denominator)


def product_solution(r, x0, t, sigma):
    return (
        free_particle_analytical_solution(r[0], x0, t, sigma)
        * free_particle_analytical_solution(r[1], x0, t, sigma)
        * free_particle_analytical_solution(r[2], x0, t, sigma)
    )


def relative_error(mra, precision, sigma, x0, t1, t, finest_scale=vp3.TimeEvolutionOperator.Adaptive):
    Pc = vp3.ScalingProjector(mra, prec=precision, dtype=complex)
    f = Pc(lambda r: product_solution(r, x0, t1, sigma))
    g = Pc(lambda r: product_solution(r, x0, t1 + t, sigma))
    U = vp3.TimeEvolutionOperator(mra, precision, t, finest_scale)
    assert U.iscomplex()
    out = U(f)
    assert isinstance(out, vp3.ComplexFunctionTree)
    return (out - g).norm() / g.norm()


def test_time_evolution_unit_box():
    mra = vp3.MultiResolutionAnalysis(vp3.BoundingBox(0), lbasis)
    precision = 1.0e-4
    assert relative_error(mra, precision, 0.02, 0.5, 0.001, 0.03, finest_scale=7) < precision / 10  # 5.1e-6


def test_time_evolution_on_box():
    """[-30, 30]^3: root scale -4, two root boxes per direction."""
    mra = vp3.MultiResolutionAnalysis(vp3.BoundingBox([-30, 30]), lbasis)
    precision = 1.0e-5
    assert relative_error(mra, precision, 50.0, 0.5, 0.0, 100.0) < precision / 10  # 8.1e-8


def test_real_input_is_promoted():
    mra = vp3.MultiResolutionAnalysis(vp3.BoundingBox([-30, 30]), lbasis)
    precision = 1.0e-5
    P = vp3.ScalingProjector(mra, prec=precision)
    f = P(vp3.GaussFunc(beta=1.0 / 50.0, alpha=1.0, position=[0.5, 0.5, 0.5]))

    U = vp3.TimeEvolutionOperator(mra, precision, 100.0)
    from_real = U(f)
    from_cplx = U(f.to_complex())
    assert isinstance(from_real, vp3.ComplexFunctionTree)
    assert (from_real - from_cplx).norm() / from_cplx.norm() == pytest.approx(0.0, abs=1e-12)
