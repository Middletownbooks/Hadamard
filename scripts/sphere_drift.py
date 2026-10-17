"""Sphere model of the partial-join drift (energy shell): four independent uniform points on spheres S^{n-1}(R_k)
(fixed PAF energy), density at 0 of the sum of their projections on j coordinates, against the Gaussian with the same
covariance; prints log2(p_gauss / p_sphere) per j.  usage: python3 scripts/sphere_drift.py N R2_A,R2_B,R2_C,R2_D J1,J2,..."""


import numpy as np, sys
from scipy.special import jv, gammaln
from scipy.integrate import quad
def logphi(r, R, n):           # log |ch.f.| of a uniform point on S^{n-1}(R) at |t| = r (signed handled separately)
    nu = n / 2 - 1; x = r * R
    return gammaln(n / 2) + nu * np.log(2 / x)
def phi(r, R, n):
    nu = n / 2 - 1; x = r * R
    if x < 1e-8: return 1.0
    return np.exp(gammaln(n / 2) + nu * np.log(2 / x)) * jv(nu, x)
def excess(j, Rs, n):
    # log2 p_gauss(0) - log2 p_sphere(0) for the j-dim projections of the sum
    var = sum(R * R / n for R in Rs)
    f = lambda r: np.prod([phi(r, R, n) for R in Rs]) * r ** (j - 1)
    rmax = 60 / min(Rs)
    I = sum(quad(f, a, b, limit=400)[0] for a, b in zip(np.linspace(0, rmax, 61)[:-1], np.linspace(0, rmax, 61)[1:]))
    g = lambda r: np.exp(-var * r * r / 2) * r ** (j - 1)
    Ig = quad(g, 0, rmax, limit=400)[0]
    return float(np.log2(Ig / I))
n = int(sys.argv[1]); Rs = [float(x) ** .5 for x in sys.argv[2].split(',')]
print('n', n, 'R^2', [round(R * R) for R in Rs])
for j in [int(x) for x in sys.argv[3].split(',')]:
    print(j, round(excess(j, Rs, n), 3))
