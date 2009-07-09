"""
Number theory module (primes, etc)
"""

from .generate import nextprime, prevprime, prime, primepi, primerange,\
randprime, Sieve, sieve
from .primetest import isprime, mr, mr_safe
from .factor_ import divisors, factorint, multiplicity, perfect_power,\
pollard_pm1, pollard_rho, primefactors, totient, trailing
from .partitions_ import npartitions
from .residue import is_primitive_root, is_quad_residue, legendre_symbol,\
n_order, totient_
from .multinomial import binomial_coefficients, binomial_coefficients_list,\
multinomial_coefficients
