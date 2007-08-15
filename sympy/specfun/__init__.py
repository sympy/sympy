"""
This module contains common special functions such
as trigonometric functions, orthogonal polynomials, the gamma function,
and so on.
"""

from factorials import factorial, binomial2, rising_factorial, \
    falling_factorial, gamma, lower_gamma, upper_gamma, \
    factorial_simplify

from orthogonal_polynomials import legendre, legendre_zero, \
    chebyshev_zero

from combinatorial import fibonacci, lucas, bernoulli, bell, harmonic

from zeta_functions import zeta, dirichlet_eta, polygamma, digamma, \
    trigamma, tetragamma
