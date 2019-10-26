"""
SymPy statistics module

Introduces a random variable type into the SymPy language.

Random variables may be declared using prebuilt functions such as
Normal, Exponential, Coin, Die, etc...  or built with functions like FiniteRV.

Queries on random expressions can be made using the functions

========================= =============================
    Expression                    Meaning
------------------------- -----------------------------
 ``P(condition)``          Probability
 ``E(expression)``         Expected value
 ``H(expression)``         Entropy
 ``variance(expression)``  Variance
 ``density(expression)``   Probability Density Function
 ``sample(expression)``    Produce a realization
 ``where(condition)``      Where the condition is true
========================= =============================

Examples
========

>>> from sympy.stats import P, E, variance, Die, Normal
>>> from sympy import Eq, simplify
>>> X, Y = Die('X', 6), Die('Y', 6) # Define two six sided dice
>>> Z = Normal('Z', 0, 1) # Declare a Normal random variable with mean 0, std 1
>>> P(X>3) # Probability X is greater than 3
1/2
>>> E(X+Y) # Expectation of the sum of two dice
7
>>> variance(X+Y) # Variance of the sum of two dice
35/6
>>> simplify(P(Z>1)) # Probability of Z being greater than 1
1/2 - erf(sqrt(2)/2)/2
"""

__all__ = []

from . import rv_interface
from .rv_interface import (
    cdf, characteristic_function, covariance, density, dependent, E, given, independent, P, pspace,
    random_symbols, sample, sample_iter, skewness, kurtosis, std, variance, where, factorial_moment,
    correlation, moment, cmoment, smoment, sampling_density, moment_generating_function, entropy, H,
    quantile
)
__all__.extend(rv_interface.__all__)

from . import frv_types
from .frv_types import (
    Bernoulli, Binomial, BetaBinomial, Coin, Die, DiscreteUniform, FiniteRV, Hypergeometric,
    Rademacher,
)
__all__.extend(frv_types.__all__)

from . import crv_types
from .crv_types import (
    ContinuousRV,
    Arcsin, Benini, Beta, BetaNoncentral, BetaPrime, Cauchy, Chi, ChiNoncentral, ChiSquared,
    Dagum, Erlang, ExGaussian, Exponential, ExponentialPower, FDistribution, FisherZ, Frechet,
    Gamma, GammaInverse, Gumbel, Gompertz, Kumaraswamy, Laplace, Logistic, LogLogistic, LogNormal,
    Maxwell, Nakagami, Normal, GaussianInverse, Pareto, QuadraticU, RaisedCosine, Rayleigh,
    ShiftedGompertz, StudentT, Trapezoidal, Triangular, Uniform, UniformSum, VonMises,
    Weibull, WignerSemicircle, Wald
)
__all__.extend(crv_types.__all__)

from . import drv_types
from .drv_types import (Geometric, Logarithmic, NegativeBinomial, Poisson, Skellam,
    YuleSimon, Zeta)
__all__.extend(drv_types.__all__)

from . import joint_rv_types
from .joint_rv_types import (
    JointRV,
    Dirichlet, GeneralizedMultivariateLogGamma, GeneralizedMultivariateLogGammaOmega,
    Multinomial, MultivariateBeta, MultivariateEwens, MultivariateT, NegativeMultinomial,
    NormalGamma
)
__all__.extend(joint_rv_types.__all__)

from . import stochastic_process_types
from .stochastic_process_types import (
    StochasticProcess,
    ContinuousTimeStochasticProcess,
    DiscreteTimeStochasticProcess,
    DiscreteMarkovChain,
    TransitionMatrixOf,
    StochasticStateSpaceOf,
    ContinuousMarkovChain,
    GeneratorMatrixOf
)
__all__.extend(stochastic_process_types.__all__)

from . import random_matrix_models
from .random_matrix_models import (
    CircularEnsemble,
    CircularUnitaryEnsemble,
    CircularOrthogonalEnsemble,
    CircularSymplecticEnsemble,
    GaussianEnsemble,
    GaussianUnitaryEnsemble,
    GaussianOrthogonalEnsemble,
    GaussianSymplecticEnsemble,
    JointEigenDistribution,
    joint_eigen_distribution,
    level_spacing_distribution
)
__all__.extend(random_matrix_models.__all__)

from . import symbolic_probability
from .symbolic_probability import Probability, Expectation, Variance, Covariance
__all__.extend(symbolic_probability.__all__)
