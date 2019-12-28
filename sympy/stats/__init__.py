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

__all__ = [
    'P', 'E', 'H', 'density', 'where', 'given', 'sample', 'cdf',
    'characteristic_function', 'pspace', 'sample_iter', 'variance', 'std',
    'skewness', 'kurtosis', 'covariance', 'dependent', 'entropy', 'independent',
    'random_symbols', 'correlation', 'factorial_moment', 'moment', 'cmoment',
    'sampling_density', 'moment_generating_function', 'smoment', 'quantile',

    'FiniteRV', 'DiscreteUniform', 'Die', 'Bernoulli', 'Coin', 'Binomial',
    'BetaBinomial', 'Hypergeometric', 'Rademacher',

    'ContinuousRV', 'Arcsin', 'Benini', 'Beta', 'BetaNoncentral', 'BetaPrime',
    'Cauchy', 'Chi', 'ChiNoncentral', 'ChiSquared', 'Dagum', 'Erlang',
    'ExGaussian', 'Exponential', 'ExponentialPower', 'FDistribution',
    'FisherZ', 'Frechet', 'Gamma', 'GammaInverse', 'Gompertz', 'Gumbel',
    'Kumaraswamy', 'Laplace', 'Logistic', 'LogLogistic', 'LogNormal',
    'Maxwell', 'Nakagami', 'Normal', 'GaussianInverse', 'Pareto',
    'PowerFunction', 'QuadraticU', 'RaisedCosine', 'Rayleigh', 'StudentT',
    'ShiftedGompertz', 'Trapezoidal', 'Triangular', 'Uniform', 'UniformSum',
    'VonMises', 'Wald', 'Weibull', 'WignerSemicircle',

    'Geometric', 'Logarithmic', 'NegativeBinomial', 'Poisson', 'Skellam',
    'YuleSimon', 'Zeta',

    'JointRV', 'Dirichlet', 'GeneralizedMultivariateLogGamma',
    'GeneralizedMultivariateLogGammaOmega', 'Multinomial', 'MultivariateBeta',
    'MultivariateEwens', 'MultivariateT', 'NegativeMultinomial',
    'NormalGamma',

    'StochasticProcess', 'DiscreteTimeStochasticProcess',
    'DiscreteMarkovChain', 'TransitionMatrixOf', 'StochasticStateSpaceOf',
    'GeneratorMatrixOf', 'ContinuousMarkovChain',

    'CircularEnsemble', 'CircularUnitaryEnsemble',
    'CircularOrthogonalEnsemble', 'CircularSymplecticEnsemble',
    'GaussianEnsemble', 'GaussianUnitaryEnsemble',
    'GaussianOrthogonalEnsemble', 'GaussianSymplecticEnsemble',
    'joint_eigen_distribution', 'JointEigenDistribution',
    'level_spacing_distribution',

    'Probability', 'Expectation', 'Variance', 'Covariance',

]
from .rv_interface import (P, E, H, density, where, given, sample, cdf,
        characteristic_function, pspace, sample_iter, variance, std, skewness,
        kurtosis, covariance, dependent, entropy, independent, random_symbols,
        correlation, factorial_moment, moment, cmoment, sampling_density,
        moment_generating_function, smoment, quantile)

from .frv_types import (FiniteRV, DiscreteUniform, Die, Bernoulli, Coin,
        Binomial, BetaBinomial, Hypergeometric, Rademacher)

from .crv_types import (ContinuousRV, Arcsin, Benini, Beta, BetaNoncentral,
        BetaPrime, Cauchy, Chi, ChiNoncentral, ChiSquared, Dagum, Erlang,
        ExGaussian, Exponential, ExponentialPower, FDistribution, FisherZ,
        Frechet, Gamma, GammaInverse, Gompertz, Gumbel, Kumaraswamy, Laplace,
        Logistic, LogLogistic, LogNormal, Maxwell, Nakagami, Normal,
        GaussianInverse, Pareto, PowerFunction, QuadraticU, RaisedCosine, Rayleigh,
        StudentT,ShiftedGompertz, Trapezoidal, Triangular, Uniform, UniformSum,
        VonMises, Wald, Weibull, WignerSemicircle)

from .drv_types import (Geometric, Logarithmic, NegativeBinomial, Poisson,
        Skellam, YuleSimon, Zeta)

from .joint_rv_types import (JointRV, Dirichlet,
        GeneralizedMultivariateLogGamma, GeneralizedMultivariateLogGammaOmega,
        Multinomial, MultivariateBeta, MultivariateEwens, MultivariateT,
        NegativeMultinomial, NormalGamma)

from .stochastic_process_types import (StochasticProcess,
        DiscreteTimeStochasticProcess, DiscreteMarkovChain,
        TransitionMatrixOf, StochasticStateSpaceOf, GeneratorMatrixOf,
        ContinuousMarkovChain)

from .random_matrix_models import (CircularEnsemble, CircularUnitaryEnsemble,
        CircularOrthogonalEnsemble, CircularSymplecticEnsemble,
        GaussianEnsemble, GaussianUnitaryEnsemble, GaussianOrthogonalEnsemble,
        GaussianSymplecticEnsemble, joint_eigen_distribution,
        JointEigenDistribution, level_spacing_distribution)

from .symbolic_probability import (Probability, Expectation, Variance,
        Covariance)
