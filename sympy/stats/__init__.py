"""
SymPy statistics module
"""

from rv import pspace, random_symbols
from rv_interface import (P, E, Density, Where, Given, CDF, Var, Std, Covar,
        Skewness, dependent, independent, Sample, sample_iter)
from frv_examples import Die, Bernoulli, Coin, FiniteRV
from crv_examples import Normal, Exponential, Gamma, Beta, Pareto, Uniform
