"""
This module contains common special functions such
as trigonometric functions, orthogonal polynomials, the gamma function,
and so on.
"""

from gamma_functions import gamma, lowergamma, uppergamma

from factorials import factorial, binomial2, rising_factorial, \
    falling_factorial, factorial_simplify

from orthogonal_polynomials import legendre, \
    chebyshev_zero

from combinatorial import fibonacci, lucas, bernoulli, bell, harmonic

from zeta_functions import zeta, dirichlet_eta, polygamma, digamma, \
    trigamma, tetragamma
