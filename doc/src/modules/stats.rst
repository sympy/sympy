=====
Stats
=====

.. automodule:: sympy.stats

Random Variable Types
^^^^^^^^^^^^^^^^^^^^^

Finite Types
------------
.. autofunction:: DiscreteUniform
.. autofunction:: Die
.. autofunction:: Bernoulli
.. autofunction:: Coin
.. autofunction:: Binomial
.. autofunction:: BetaBinomial
.. autofunction:: Hypergeometric
.. autofunction:: FiniteRV
.. autofunction:: Rademacher

Discrete Types
--------------
.. autofunction:: Geometric
.. autofunction:: Hermite
.. autofunction:: Poisson
.. autofunction:: Logarithmic
.. autofunction:: NegativeBinomial
.. autofunction:: Skellam
.. autofunction:: YuleSimon
.. autofunction:: Zeta

Continuous Types
----------------

.. autofunction:: Arcsin
.. autofunction:: Benini
.. autofunction:: Beta
.. autofunction:: BetaNoncentral
.. autofunction:: BetaPrime
.. autofunction:: BoundedPareto
.. autofunction:: Cauchy
.. autofunction:: Chi
.. autofunction:: ChiNoncentral
.. autofunction:: ChiSquared
.. autofunction:: Dagum
.. autofunction:: Davis
.. autofunction:: Erlang
.. autofunction:: ExGaussian
.. autofunction:: Exponential
.. autofunction:: FDistribution
.. autofunction:: FisherZ
.. autofunction:: Frechet
.. autofunction:: Gamma
.. autofunction:: GammaInverse
.. autofunction:: Gompertz
.. autofunction:: Gumbel
.. autofunction:: Kumaraswamy
.. autofunction:: Laplace
.. autofunction:: Levy
.. autofunction:: Logistic
.. autofunction:: LogLogistic
.. autofunction:: LogNormal
.. autofunction:: Lomax
.. autofunction:: Maxwell
.. autofunction:: Moyal
.. autofunction:: Nakagami
.. autofunction:: Normal
.. autofunction:: Pareto
.. autofunction:: PowerFunction
.. autofunction:: QuadraticU
.. autofunction:: RaisedCosine
.. autofunction:: Rayleigh
.. autofunction:: Reciprocal
.. autofunction:: StudentT
.. autofunction:: ShiftedGompertz
.. autofunction:: Trapezoidal
.. autofunction:: Triangular
.. autofunction:: Uniform
.. autofunction:: UniformSum
.. autofunction:: VonMises
.. autofunction:: Wald
.. autofunction:: Weibull
.. autofunction:: WignerSemicircle
.. autofunction:: ContinuousRV

Joint Types
-----------
.. autofunction:: JointRV
.. autofunction:: marginal_distribution
.. autofunction:: MultivariateNormal
.. autofunction:: MultivariateLaplace
.. autofunction:: GeneralizedMultivariateLogGamma
.. autofunction:: GeneralizedMultivariateLogGammaOmega
.. autofunction:: Multinomial
.. autofunction:: MultivariateBeta
.. autofunction:: MultivariateEwens
.. autofunction:: MultivariateT
.. autofunction:: NegativeMultinomial
.. autofunction:: NormalGamma

.. _sympy-stats-stochastic-processes:

Stochastic Processes
--------------------
.. autoclass:: DiscreteMarkovChain
   :members:
.. autoclass:: ContinuousMarkovChain
   :members:
.. autoclass:: BernoulliProcess
   :members:
.. autoclass:: PoissonProcess
   :members:
.. autoclass:: WienerProcess
   :members:
.. autoclass:: GammaProcess
   :members:

Matrix Distributions
--------------------
.. autofunction:: MatrixGamma
.. autofunction:: Wishart
.. autofunction:: MatrixNormal

Compound Distribution
---------------------
.. autoclass:: sympy.stats.compound_rv.CompoundDistribution
   :members:

Interface
^^^^^^^^^

.. autofunction:: P
.. autoclass:: Probability
   :members:
.. autofunction:: E
.. autoclass:: Expectation
   :members:
.. autofunction:: density
.. autofunction:: entropy
.. autofunction:: given
.. autofunction:: where
.. autofunction:: variance
.. autoclass:: Variance
   :members:
.. autofunction:: covariance
.. autoclass:: Covariance
   :members:
.. autofunction:: coskewness
.. autofunction:: median
.. autofunction:: std
.. autofunction:: quantile
.. autofunction:: sample
.. autofunction:: sample_iter
.. autofunction:: factorial_moment
.. autofunction:: kurtosis
.. autofunction:: skewness
.. autofunction:: correlation
.. autofunction:: sympy.stats.rv.sampling_density
.. autofunction:: sympy.stats.rv.sampling_P
.. autofunction:: sympy.stats.rv.sampling_E
.. autoclass:: Moment
   :members:
.. autofunction:: moment
.. autoclass:: CentralMoment
   :members:
.. autofunction:: cmoment
.. autoclass:: ExpectationMatrix
   :members:
.. autoclass:: VarianceMatrix
   :members:
.. autoclass:: CrossCovarianceMatrix
   :members:

Mechanics
^^^^^^^^^
.. module:: sympy.stats.rv

SymPy Stats employs a relatively complex class hierarchy.

``RandomDomain``\s are a mapping of variables to possible values.
For example, we might say that the symbol ``Symbol('x')`` can
take on the values `\{1,2,3,4,5,6\}`.

.. class:: RandomDomain

A ``PSpace``, or Probability Space, combines a ``RandomDomain``
with a density to provide probabilistic information. For example
the above domain could be enhanced by a finite density
``{1:1/6, 2:1/6, 3:1/6, 4:1/6, 5:1/6, 6:1/6}`` to
fully define the roll of a fair die named ``x``.

.. class:: PSpace

A RandomSymbol represents the PSpace's symbol 'x' inside of SymPy expressions.

.. class:: RandomSymbol

The RandomDomain and PSpace classes are almost never directly instantiated.
Instead they are subclassed for a variety of situations.

RandomDomains and PSpaces must be sufficiently general to represent domains and
spaces of several variables with arbitrarily complex densities. This generality
is often unnecessary. Instead we often build SingleDomains and SinglePSpaces to
represent single, univariate events and processes such as a single die or a
single normal variable.

.. class:: SinglePSpace
.. class:: SingleDomain


Another common case is to collect together a set of such univariate random
variables. A collection of independent SinglePSpaces or SingleDomains can be
brought together to form a ProductDomain or ProductPSpace. These objects would
be useful in representing three dice rolled together for example.

.. class:: ProductDomain

.. class:: ProductPSpace

The Conditional adjective is added whenever we add a global condition to a
RandomDomain or PSpace. A common example would be three independent dice where
we know their sum to be greater than 12.

.. class:: ConditionalDomain

We specialize further into Finite and Continuous versions of these classes to
represent finite (such as dice) and continuous (such as normals) random
variables.

.. module:: sympy.stats.frv
.. class:: FiniteDomain
.. class:: FinitePSpace

.. module:: sympy.stats.crv
.. class:: ContinuousDomain
.. class:: ContinuousPSpace

Additionally there are a few specialized classes that implement certain common
random variable types. There is for example a DiePSpace that implements
SingleFinitePSpace and a NormalPSpace that implements SingleContinuousPSpace.

.. module:: sympy.stats.frv_types
.. class:: DiePSpace

.. module:: sympy.stats.crv_types
.. class:: NormalPSpace

RandomVariables can be extracted from these objects using the PSpace.values
method.

As previously mentioned SymPy Stats employs a relatively complex class
structure. Inheritance is widely used in the implementation of end-level
classes. This tactic was chosen to balance between the need to allow SymPy to
represent arbitrarily defined random variables and optimizing for common cases.
This complicates the code but is structured to only be important to those
working on extending SymPy Stats to other random variable types.

Users will not use this class structure. Instead these mechanics are exposed
through variable creation functions Die, Coin, FiniteRV, Normal, Exponential,
etc.... These build the appropriate SinglePSpaces and return the corresponding
RandomVariable. Conditional and Product spaces are formed in the natural
construction of SymPy expressions and the use of interface functions E, Given,
Density, etc....


.. function:: sympy.stats.Die
.. function:: sympy.stats.Normal

There are some additional functions that may be useful. They are largely used
internally.


.. autofunction:: sympy.stats.rv.random_symbols
.. autofunction:: sympy.stats.rv.pspace
.. autofunction:: sympy.stats.rv.rs_swap
