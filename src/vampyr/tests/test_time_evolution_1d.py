from vampyr import vampyr1d as vp1
import numpy as np
import pytest

#Analytical solution parameters
x0 = 0.5
sigma = 0.001

#Time step
time = 0.002

# Set the precision and make the MRA
order = 5
precision = 1.0e-6
finest_scale = 8


from vampyr import InterpolatingBasis, LegendreBasis
lbasis = LegendreBasis(order=order)
box = vp1.BoundingBox(0)
mra = vp1.MultiResolutionAnalysis(box, lbasis) 

# Make the scaling projectors
P = vp1.ScalingProjector(mra, prec = precision)
Pc = vp1.ComplexScalingProjector(mra, prec = precision)

# Make the Schrodinger free-particle semigroup operator
semigroup_operator = vp1.TimeEvolutionOperator(mra, precision, time, finest_scale)

def free_particle_analytical_solution(x, x0, t, sigma):
    denominator = 4.0j * t + sigma
    sqrt_denom = np.sqrt(denominator)
    exponent = -((x - x0) ** 2) / denominator
    return np.sqrt(sigma) / sqrt_denom * np.exp(exponent)

def f(x):
    return np.exp( - (x[0] - x0)**2 / sigma )
f = P(f)

def g(x):
    return free_particle_analytical_solution(x[0], x0, time, sigma)
g = Pc(g)

# The input is real, the result complex
out = semigroup_operator(f)

difference = out - g                    #1.4e-15

#print(difference.squaredNorm())

epsilon = 5e-14

def test_time_evolution():
    assert semigroup_operator.iscomplex()
    assert isinstance(out, vp1.ComplexFunctionTree)
    assert difference.squaredNorm() == pytest.approx(0.0, abs = epsilon)

def test_time_evolution_parts():
    assert (out.real() - g.real()).squaredNorm() == pytest.approx(0.0, abs = epsilon)
    assert (out.imag() - g.imag()).squaredNorm() == pytest.approx(0.0, abs = epsilon)

def test_time_evolution_two_steps():
    # After the first step the packet has reached the box edges, hence the
    # looser bound
    out2 = semigroup_operator(out)
    def g2(x):
        return free_particle_analytical_solution(x[0], x0, 2 * time, sigma)
    assert (out2 - Pc(g2)).squaredNorm() == pytest.approx(0.0, abs = 1e-6)
    # The norm is conserved up to what has left the box
    assert out.norm() == pytest.approx(f.norm(), rel = 1e-4)

def test_real_input_is_promoted():
    # MRCPP aborts on a real tree, so the binding promotes it first
    from_complex = semigroup_operator(f.to_complex())
    assert (out - from_complex).norm() / from_complex.norm() == pytest.approx(0.0, abs = 1e-10)
    assert out.imag().squaredNorm() > 1e-12

def test_time_evolution_adaptive():
    # The default build is adaptive
    adaptive = vp1.TimeEvolutionOperator(mra, precision, time)
    assert (adaptive(f) - g).squaredNorm() == pytest.approx(0.0, abs = 1e-12)

def test_time_evolution_on_box():
    """[-30, 30]: root scale -4, two root boxes, scaling factor 1.875."""
    mra_box = vp1.MultiResolutionAnalysis(vp1.BoundingBox([-30, 30]), lbasis)
    Pc_box = vp1.ComplexScalingProjector(mra_box, prec = precision)

    def relative_error(sigma, t, finest_scale):
        U = vp1.TimeEvolutionOperator(mra_box, precision, t, finest_scale)
        f = Pc_box(lambda x: complex(np.exp(-(x[0] - x0)**2 / sigma)))
        g = Pc_box(lambda x: free_particle_analytical_solution(x[0], x0, t, sigma))
        return (U(f) - g).norm() / g.norm()

    assert relative_error(5.0, 10.0, vp1.TimeEvolutionOperator.Adaptive) < precision   #1.9e-7
    assert relative_error(50.0, 0.1, 2) < precision                                     #4.4e-8

def test_bad_arguments_raise():
    U = vp1.TimeEvolutionOperator
    with pytest.raises(ValueError):
        U(vp1.MultiResolutionAnalysis(box, InterpolatingBasis(order=order)), precision, time)
    with pytest.raises(ValueError):
        U(vp1.MultiResolutionAnalysis(vp1.BoundingBox(scaling=[1.0], pbc=True), lbasis), precision, time)
    for bad_time in (0.0, -time):
        with pytest.raises(ValueError):
            U(mra, precision, bad_time)
    for bad_Jpower in (0, 120):
        with pytest.raises(ValueError):
            U(mra, precision, time, max_Jpower = bad_Jpower)
    # -1 is coarser than the root scale 0 of the unit box, but finer than the
    # root scale -4 of [-30, 30]
    with pytest.raises(ValueError):
        U(mra, precision, time, finest_scale = -1)
    U(vp1.MultiResolutionAnalysis(vp1.BoundingBox([-30, 30]), lbasis), precision, 0.1, finest_scale = -1)

def test_advanced_apply_needs_complex_trees():
    with pytest.raises(TypeError):
        vp1.advanced.apply(precision, vp1.FunctionTree(mra), semigroup_operator, f)

def test_max_Jpower_reaches_the_adaptive_path():
    full = vp1.TimeEvolutionOperator(mra, precision, time, max_Jpower = 30)
    short = vp1.TimeEvolutionOperator(mra, precision, time, max_Jpower = 2)
    assert (full(f) - short(f)).norm() > 1e-10

def test_legacy_imaginary_flag_is_rejected():
    # Without these overloads the flag would bind as finest_scale or max_Jpower
    with pytest.raises(TypeError):
        vp1.TimeEvolutionOperator(mra, precision, time, True)
    with pytest.raises(TypeError):
        vp1.TimeEvolutionOperator(mra, precision, time, finest_scale, False)
