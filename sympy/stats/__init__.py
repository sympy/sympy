"""
SymPy statistics module

Introduces a random variable type into the SymPy language.

Random variables may be declared using prebuilt functions such as
Normal, Exponential, Coin, Die, etc...  or built with functions like FiniteRV.

Queries on random expressions can be made using the functions

P(condition) -- Probability
E(expression) -- Expectation value
Var(expression) -- Variance
Density(expression) -- Probability Density Function (PDF)
Sample(expression) -- Produce a realization
Where(condition) -- Domain of where the condition might be true

>>> from sympy.stats import P, E, Var, Die, Normal
>>> from sympy import Eq
>>> X, Y = Die(6), Die(6) # Define two six sided dice
>>> Y = Normal(0, 1) # Declare a Normal random variable with mean 0, std 1
>>> P(X>3) # Probability X is greater than 3
1/2
>>> E(X+Y) # Expectation of the sum of two dice
7
>>> Var(X+Y) # Variance of the sum of two dice
35/6
>>> P(Z>1) # Probability of Z being greater than 1
-erf(sqrt(2)/2)/2 + 1/2
"""

from rv import pspace, random_symbols
from rv_interface import (P, E, Density, Where, Given, CDF, Var, Std, Covar,
        Skewness, dependent, independent, Sample, sample_iter)
from frv_examples import Die, Bernoulli, Coin, FiniteRV
from crv_examples import Normal, Exponential, Gamma, Beta, Pareto, Uniform
