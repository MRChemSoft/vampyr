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

# Make the scaling projector
P = vp1.ScalingProjector(mra, prec = precision)

# Make the Schrodinger free-particle semigroup operator
imaginary = False
real_semigroup_operator = vp1.TimeEvolutionOperator(mra, precision, time, 
finest_scale, imaginary)
imaginary = True
imag_semigroup_operator = vp1.TimeEvolutionOperator(mra, precision, time, 
finest_scale, imaginary)

def free_particle_analytical_solution(x, x0, t, sigma):
    denominator = 4.0j * t + sigma
    sqrt_denom = np.sqrt(denominator)
    exponent = -((x - x0) ** 2) / denominator
    return np.sqrt(sigma) / sqrt_denom * np.exp(exponent)

def f(x):
    return np.exp( - (x[0] - x0)**2 / sigma )
f = P(f)

def Re_g(x):
    return free_particle_analytical_solution(x[0], x0, time, sigma).real
Re_g = P(Re_g)

def Im_g(x):
    return free_particle_analytical_solution(x[0], x0, time, sigma).imag
Im_g = P(Im_g)

Re_out = real_semigroup_operator(f)
Im_out = imag_semigroup_operator(f)

Re_difference = Re_out - Re_g           #2.1e-14
Im_difference = Im_out - Im_g           #2.2e-14

#print(Re_difference.squaredNorm())
#print(Im_difference.squaredNorm())

epsilon = 2.5e-14

assert Re_difference.squaredNorm() == pytest.approx(0.0, abs = epsilon)
assert Im_difference.squaredNorm() == pytest.approx(0.0, abs = epsilon)