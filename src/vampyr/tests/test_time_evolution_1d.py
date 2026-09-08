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


from vampyr import LegendreBasis
lbasis = LegendreBasis(order=order)
box = vp1.BoundingBox(0)
mra = vp1.MultiResolutionAnalysis(box, lbasis) 

# Make the scaling projectors
P = vp1.ScalingProjector(mra, prec = precision)
Pc = vp1.ScalingProjector(mra, prec = precision, dtype = complex)

# The kernel is complex and carries the whole cos + 1j*sin coefficient, so this
# is one operator, not a real/imaginary pair.
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

# The real input is promoted on the way in, so the result is complex.
out = semigroup_operator(f)

difference = out - g                    #4.3e-14

#print(difference.squaredNorm())

epsilon = 5e-14

def test_time_evolution():
    assert semigroup_operator.iscomplex()
    assert isinstance(out, vp1.ComplexFunctionTree)
    assert difference.squaredNorm() == pytest.approx(0.0, abs = epsilon)

def test_time_evolution_parts():
    # What the two-operator API returned separately is now the real and
    # imaginary part of one tree.
    assert (out.real() - g.real()).squaredNorm() == pytest.approx(0.0, abs = epsilon)
    assert (out.imag() - g.imag()).squaredNorm() == pytest.approx(0.0, abs = epsilon)

def test_time_evolution_adaptive():
    # Adaptive is the default build now.
    adaptive = vp1.TimeEvolutionOperator(mra, precision, time)
    assert (adaptive(f) - g).squaredNorm() == pytest.approx(0.0, abs = 1e-12)
