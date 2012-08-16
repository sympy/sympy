"""
Continuous Random Variables - Prebuilt variables

Contains
========
Arcsin
Benini
Beta
BetaPrime
Cauchy
Chi
Dagum
Exponential
Gamma
Laplace
Logistic
LogNormal
Maxwell
Nakagami
Normal
Pareto
Rayleigh
StudentT
Triangular
Uniform
UniformSum
Weibull
WignerSemicircle
"""

from sympy import (exp, log, sqrt, pi, S, Dummy, Interval, S, sympify, gamma,
                   Piecewise, And, Eq, binomial, factorial, Sum, floor, Abs,
                   Symbol, log)
from sympy import beta as beta_fn
from crv import SingleContinuousPSpace
from sympy.core.decorators import _sympifyit
import random

oo = S.Infinity

__all__ = ['ContinuousRV',
'Arcsin',
'Benini',
'Beta',
'BetaPrime',
'Cauchy',
'Chi',
'Dagum',
'Exponential',
'Gamma',
'Laplace',
'Logistic',
'LogNormal',
'Maxwell',
'Nakagami',
'Normal',
'Pareto',
'Rayleigh',
'StudentT',
'Triangular',
'Uniform',
'UniformSum',
'Weibull',
'WignerSemicircle'
]

def _value_check(condition, message):
    """
    Check a condition on input value.

    Raises ValueError with message if condition is not True
    """
    if condition is not True:
        raise ValueError(message)


def ContinuousRV(symbol, density, set=Interval(-oo,oo)):
    """
    Create a Continuous Random Variable given the following:

    -- a symbol
    -- a probability density function
    -- set on which the pdf is valid (defaults to entire real line)

    Returns a RandomSymbol.

    Many common continuous random variable types are already implemented.
    This function should be necessary only very rarely.

    Examples
    ========

    >>> from sympy import Symbol, sqrt, exp, pi
    >>> from sympy.stats import ContinuousRV, P, E

    >>> x = Symbol("x")

    >>> pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)) # Normal distribution
    >>> X = ContinuousRV(x, pdf)

    >>> E(X)
    0
    >>> P(X>0)
    1/2
    """
    return SingleContinuousPSpace(symbol, density, set).value

########################################
# Continuous Probability Distributions #
########################################

#-------------------------------------------------------------------------------
# Arcsin distribution ----------------------------------------------------------

class ArcsinPSpace(SingleContinuousPSpace):
    def __new__(cls, name, a, b):
        a, b = sympify(a), sympify(b)
        x = Symbol(name)
        pdf = 1/(pi*sqrt((x-a)*(b-x)))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(a, b))
        return obj

def Arcsin(name, a=0, b=1):
    r"""
    Create a Continuous Random Variable with an arcsin distribution.

    The density of the arcsin distribution is given by

    .. math::
        f(x) := \frac{1}{\pi\sqrt{(x-a)(b-x)}}

    with :math:`x \in [a,b]`. It must hold that :math:`-\infty < a < b < \infty`.

    Parameters
    ==========

    a : Real number, the left interval boundary
    b : Real number, the right interval boundary

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Arcsin, density
    >>> from sympy import Symbol, simplify

    >>> a = Symbol("a", real=True)
    >>> b = Symbol("b", real=True)

    >>> X = Arcsin("x", a, b)

    >>> density(X)
    Lambda(_x, 1/(pi*sqrt((-_x + b)*(_x - a))))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Arcsine_distribution
    """

    return ArcsinPSpace(name, a, b).value

#-------------------------------------------------------------------------------
# Benini distribution ----------------------------------------------------------

class BeniniPSpace(SingleContinuousPSpace):
    def __new__(cls, name, alpha, beta, sigma):
        alpha, beta, sigma = sympify(alpha), sympify(beta), sympify(sigma)
        x = Symbol(name)
        pdf = (exp(-alpha*log(x/sigma)-beta*log(x/sigma)**2)
               *(alpha/x+2*beta*log(x/sigma)/x))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf,
                                             set = Interval(sigma, oo))
        return obj

def Benini(name, alpha, beta, sigma):
    r"""
    Create a Continuous Random Variable with a Benini distribution.

    The density of the Benini distribution is given by

    .. math::
        f(x) := e^{-\alpha\log{\frac{x}{\sigma}}
                -\beta\log\left[{\frac{x}{\sigma}}\right]^2}
                \left(\frac{\alpha}{x}+\frac{2\beta\log{\frac{x}{\sigma}}}{x}\right)

    Parameters
    ==========

    alpha : Real number, `alpha` > 0 a shape
    beta : Real number, `beta` > 0 a shape
    sigma : Real number, `sigma` > 0 a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Benini, density
    >>> from sympy import Symbol, simplify, pprint

    >>> alpha = Symbol("alpha", positive=True)
    >>> beta = Symbol("beta", positive=True)
    >>> sigma = Symbol("sigma", positive=True)

    >>> X = Benini("x", alpha, beta, sigma)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /                                                             2       \
          |   /                  /  x  \\             /  x  \            /  x  \|
          |   |        2*beta*log|-----||  - alpha*log|-----| - beta*log |-----||
          |   |alpha             \sigma/|             \sigma/            \sigma/|
    Lambda|x, |----- + -----------------|*e                                     |
          \   \  x             x        /                                       /

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Benini_distribution
    """

    return BeniniPSpace(name, alpha, beta, sigma).value

#-------------------------------------------------------------------------------
# Beta distribution ------------------------------------------------------------

class BetaPSpace(SingleContinuousPSpace):
    def __new__(cls, name, alpha, beta):
        alpha, beta = sympify(alpha), sympify(beta)

        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(beta > 0, "Beta must be positive")

        x = Symbol(name)
        pdf = x**(alpha-1) * (1-x)**(beta-1) / beta_fn(alpha, beta)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, 1))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.betavariate(self.alpha, self.beta)}

def Beta(name, alpha, beta):
    r"""
    Create a Continuous Random Variable with a Beta distribution.

    The density of the Beta distribution is given by

    .. math::
        f(x) := \frac{x^{\alpha-1}(1-x)^{\beta-1}} {\mathrm{B}(\alpha,\beta)}

    with :math:`x \in [0,1]`.

    Parameters
    ==========

    alpha : Real number, `alpha` > 0 a shape
    beta : Real number, `beta` > 0 a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Beta, density, E, variance
    >>> from sympy import Symbol, simplify, pprint

    >>> alpha = Symbol("alpha", positive=True)
    >>> beta = Symbol("beta", positive=True)

    >>> X = Beta("x", alpha, beta)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /    alpha - 1         beta - 1                    \
          |   x         *(-x + 1)        *gamma(alpha + beta)|
    Lambda|x, -----------------------------------------------|
          \               gamma(alpha)*gamma(beta)           /

    >>> simplify(E(X, meijerg=True))
    alpha/(alpha + beta)

    >>> simplify(variance(X, meijerg=True))
    alpha*beta/((alpha + beta)**2*(alpha + beta + 1))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Beta_distribution
    [2] http://mathworld.wolfram.com/BetaDistribution.html
    """

    return BetaPSpace(name, alpha, beta).value

#-------------------------------------------------------------------------------
# Beta prime distribution ------------------------------------------------------

class BetaPrimePSpace(SingleContinuousPSpace):
    def __new__(cls, name, alpha, beta):
        alpha, beta = sympify(alpha), sympify(beta)
        x = Symbol(name)
        pdf = x**(alpha-1)*(1+x)**(-alpha-beta)/beta_fn(alpha, beta)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def BetaPrime(name, alpha, beta):
    r"""
    Create a continuous random variable with a Beta prime distribution.

    The density of the Beta prime distribution is given by

    .. math::
        f(x) := \frac{x^{\alpha-1} (1+x)^{-\alpha -\beta}}{B(\alpha,\beta)}

    with :math:`x > 0`.

    Parameters
    ==========

    alpha : Real number, `alpha` > 0 a shape
    beta : Real number, `beta` > 0 a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import BetaPrime, density
    >>> from sympy import Symbol, pprint

    >>> alpha = Symbol("alpha", positive=True)
    >>> beta = Symbol("beta", positive=True)

    >>> X = BetaPrime("x", alpha, beta)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /    alpha - 1        -alpha - beta                    \
          |   x         *(x + 1)             *gamma(alpha + beta)|
    Lambda|x, ---------------------------------------------------|
          \                 gamma(alpha)*gamma(beta)             /

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Beta_prime_distribution
    [2] http://mathworld.wolfram.com/BetaPrimeDistribution.html
    """

    return BetaPrimePSpace(name, alpha, beta).value

#-------------------------------------------------------------------------------
# Cauchy distribution ----------------------------------------------------------

class CauchyPSpace(SingleContinuousPSpace):
    def __new__(cls, name, x0, gamma):
        x0, gamma = sympify(x0), sympify(gamma)
        x = Symbol(name)
        pdf = 1/(pi*gamma*(1+((x-x0)/gamma)**2))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Cauchy(name, x0, gamma):
    r"""
    Create a continuous random variable with a Cauchy distribution.

    The density of the Cauchy distribution is given by

    .. math::
        f(x) := \frac{1}{\pi} \arctan\left(\frac{x-x_0}{\gamma}\right)
                +\frac{1}{2}

    Parameters
    ==========

    x0 : Real number, the location
    gamma : Real number, `gamma` > 0 the scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Cauchy, density
    >>> from sympy import Symbol

    >>> x0 = Symbol("x0")
    >>> gamma = Symbol("gamma", positive=True)

    >>> X = Cauchy("x", x0, gamma)

    >>> density(X)
    Lambda(_x, 1/(pi*gamma*(1 + (_x - x0)**2/gamma**2)))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Cauchy_distribution
    [2] http://mathworld.wolfram.com/CauchyDistribution.html
    """

    return CauchyPSpace(name, x0, gamma).value

#-------------------------------------------------------------------------------
# Chi distribution -------------------------------------------------------------

class ChiPSpace(SingleContinuousPSpace):
    def __new__(cls, name, k):
        k = sympify(k)
        x = Symbol(name)
        pdf = 2**(1-k/2)*x**(k-1)*exp(-x**2/2)/gamma(k/2)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Chi(name, k):
    r"""
    Create a continuous random variable with a Chi distribution.

    The density of the Chi distribution is given by

    .. math::
        f(x) := \frac{2^{1-k/2}x^{k-1}e^{-x^2/2}}{\Gamma(k/2)}

    with :math:`x \geq 0`.

    Parameters
    ==========

    k : Integer, `k` > 0 the number of degrees of freedom

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Chi, density, E, std
    >>> from sympy import Symbol, simplify

    >>> k = Symbol("k", integer=True)

    >>> X = Chi("x", k)

    >>> density(X)
    Lambda(_x, 2**(-k/2 + 1)*_x**(k - 1)*exp(-_x**2/2)/gamma(k/2))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Chi_distribution
    [2] http://mathworld.wolfram.com/ChiDistribution.html
    """

    return ChiPSpace(name, k).value

#-------------------------------------------------------------------------------
# Dagum distribution -----------------------------------------------------------

class DagumPSpace(SingleContinuousPSpace):
    def __new__(cls, name, p, a, b):
        p, a, b = sympify(p), sympify(a), sympify(b)
        x = Symbol(name)
        pdf = a*p/x*((x/b)**(a*p)/(((x/b)**a+1)**(p+1)))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Dagum(name, p, a, b):
    r"""
    Create a continuous random variable with a Dagum distribution.

    The density of the Dagum distribution is given by

    .. math::
        f(x) := \frac{a p}{x} \left( \frac{(\tfrac{x}{b})^{a p}}
                {\left((\tfrac{x}{b})^a + 1 \right)^{p+1}} \right)

    with :math:`x > 0`.

    Parameters
    ==========

    p : Real number, `p` > 0 a shape
    a : Real number, `a` > 0 a shape
    b : Real number, `b` > 0 a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Dagum, density
    >>> from sympy import Symbol, simplify

    >>> p = Symbol("p", positive=True)
    >>> b = Symbol("b", positive=True)
    >>> a = Symbol("a", positive=True)

    >>> X = Dagum("x", p, a, b)

    >>> density(X)
    Lambda(_x, a*p*(_x/b)**(a*p)*((_x/b)**a + 1)**(-p - 1)/_x)

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Dagum_distribution
    """

    return DagumPSpace(name, p, a, b).value

#-------------------------------------------------------------------------------
# Exponential distribution -----------------------------------------------------

class ExponentialPSpace(SingleContinuousPSpace):
    def __new__(cls, name, rate):
        rate = sympify(rate)

        _value_check(rate > 0, "Rate must be positive.")

        x = Symbol(name)
        pdf = rate * exp(-rate*x)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.rate = rate
        return obj

    def sample(self):
        return {self.value: random.expovariate(self.rate)}

def Exponential(name, rate):
    r"""
    Create a continuous random variable with an Exponential distribution.

    The density of the exponential distribution is given by

    .. math::
        f(x) := \lambda \exp(-\lambda x)

    with :math:`x > 0`.

    Parameters
    ==========

    rate : Real number, `rate` > 0 the rate or inverse scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Exponential, density, cdf, E
    >>> from sympy.stats import variance, std, skewness
    >>> from sympy import Symbol

    >>> l = Symbol("lambda", positive=True)

    >>> X = Exponential("x", l)

    >>> density(X)
    Lambda(_x, lambda*exp(-_x*lambda))

    >>> cdf(X)
    Lambda(_z, Piecewise((1 - exp(-_z*lambda), _z >= 0), (0, True)))

    >>> E(X)
    1/lambda

    >>> variance(X)
    lambda**(-2)

    >>> skewness(X)
    2

    >>> X = Exponential('x', 10)

    >>> density(X)
    Lambda(_x, 10*exp(-10*_x))

    >>> E(X)
    1/10

    >>> std(X)
    1/10

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Exponential_distribution
    [2] http://mathworld.wolfram.com/ExponentialDistribution.html
    """

    return ExponentialPSpace(name, rate).value

#-------------------------------------------------------------------------------
# Gamma distribution -----------------------------------------------------------

class GammaPSpace(SingleContinuousPSpace):
    def __new__(cls, name, k, theta):
        k, theta = sympify(k), sympify(theta)

        _value_check(k > 0, "k must be positive")
        _value_check(theta > 0, "Theta must be positive")

        x = Symbol(name)
        pdf = x**(k-1) * exp(-x/theta) / (gamma(k)*theta**k)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.k = k
        obj.theta = theta
        return obj

    def sample(self):
        return {self.value: random.gammavariate(self.k, self.theta)}

def Gamma(name, k, theta):
    r"""
    Create a continuous random variable with a Gamma distribution.

    The density of the Gamma distribution is given by

    .. math::
        f(x) := \frac{1}{\Gamma(k) \theta^k} x^{k - 1} e^{-\frac{x}{\theta}}

    with :math:`x \in [0,1]`.

    Parameters
    ==========

    k : Real number, `k` > 0 a shape
    theta : Real number, `theta` > 0 a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Gamma, density, cdf, E, variance
    >>> from sympy import Symbol, pprint

    >>> k = Symbol("k", positive=True)
    >>> theta = Symbol("theta", positive=True)

    >>> X = Gamma("x", k, theta)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /                     -x \
          |                   -----|
          |    k - 1      -k  theta|
          |   x     *theta  *e     |
    Lambda|x, ---------------------|
          \          gamma(k)      /

    >>> C = cdf(X, meijerg=True)
    >>> pprint(C, use_unicode=False)
          /   /                                   /     z  \            \
          |   |                       k*lowergamma|k, -----|            |
          |   |  k*lowergamma(k, 0)               \   theta/            |
    Lambda|z, <- ------------------ + ----------------------  for z >= 0|
          |   |     gamma(k + 1)           gamma(k + 1)                 |
          |   |                                                         |
          \   \                      0                        otherwise /

    >>> E(X)
    theta*gamma(k + 1)/gamma(k)

    >>> V = variance(X)
    >>> pprint(V, use_unicode=False)
           2      2                     -k      k + 1
      theta *gamma (k + 1)   theta*theta  *theta     *gamma(k + 2)
    - -------------------- + -------------------------------------
                2                           gamma(k)
           gamma (k)

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Gamma_distribution
    [2] http://mathworld.wolfram.com/GammaDistribution.html
    """

    return GammaPSpace(name, k, theta).value

#-------------------------------------------------------------------------------
# Laplace distribution ---------------------------------------------------------

class LaplacePSpace(SingleContinuousPSpace):
    def __new__(cls, name, mu, b):
        mu, b = sympify(mu), sympify(b)
        x = Symbol(name)
        pdf = 1/(2*b)*exp(-Abs(x-mu)/b)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Laplace(name, mu, b):
    r"""
    Create a continuous random variable with a Laplace distribution.

    The density of the Laplace distribution is given by

    .. math::
        f(x) := \frac{1}{2 b} \exp \left(-\frac{|x-\mu|}b \right)

    Parameters
    ==========

    mu : Real number, the location
    b : Real number, `b` > 0 a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Laplace, density
    >>> from sympy import Symbol

    >>> mu = Symbol("mu")
    >>> b = Symbol("b", positive=True)

    >>> X = Laplace("x", mu, b)

    >>> density(X)
    Lambda(_x, exp(-Abs(_x - mu)/b)/(2*b))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Laplace_distribution
    [2] http://mathworld.wolfram.com/LaplaceDistribution.html
    """

    return LaplacePSpace(name, mu, b).value

#-------------------------------------------------------------------------------
# Logistic distribution --------------------------------------------------------

class LogisticPSpace(SingleContinuousPSpace):
    def __new__(cls, name, mu, s):
        mu, s = sympify(mu), sympify(s)

        x = Symbol(name)
        pdf = exp(-(x-mu)/s)/(s*(1+exp(-(x-mu)/s))**2)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Logistic(name, mu, s):
    r"""
    Create a continuous random variable with a logistic distribution.

    The density of the logistic distribution is given by

    .. math::
        f(x) := \frac{e^{-(x-\mu)/s}} {s\left(1+e^{-(x-\mu)/s}\right)^2}

    Parameters
    ==========

    mu : Real number, the location
    s : Real number, `s` > 0 a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Logistic, density
    >>> from sympy import Symbol

    >>> mu = Symbol("mu", real=True)
    >>> s = Symbol("s", positive=True)

    >>> X = Logistic("x", mu, s)

    >>> density(X)
    Lambda(_x, exp((-_x + mu)/s)/(s*(exp((-_x + mu)/s) + 1)**2))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Logistic_distribution
    [2] http://mathworld.wolfram.com/LogisticDistribution.html
    """

    return LogisticPSpace(name, mu, s).value

#-------------------------------------------------------------------------------
# Log Normal distribution ------------------------------------------------------

class LogNormalPSpace(SingleContinuousPSpace):
    def __new__(cls, name, mean, std):
        mean, std = sympify(mean), sympify(std)

        x = Symbol(name)
        pdf = exp(-(log(x)-mean)**2 / (2*std**2)) / (x*sqrt(2*pi)*std)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.mean = mean
        obj.std = std
        return obj

    def sample(self):
        return {self.value: random.lognormvariate(self.mean, self.std)}

def LogNormal(name, mean, std):
    r"""
    Create a continuous random variable with a log-normal distribution.

    The density of the log-normal distribution is given by

    .. math::
        f(x) := \frac{1}{x\sqrt{2\pi\sigma^2}}
                e^{-\frac{\left(\ln x-\mu\right)^2}{2\sigma^2}}

    with :math:`x \geq 0`.

    Parameters
    ==========

    mu : Real number, the log-scale
    sigma : Real number, :math:`\sigma^2 > 0` a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import LogNormal, density
    >>> from sympy import Symbol, simplify, pprint

    >>> mu = Symbol("mu", real=True)
    >>> sigma = Symbol("sigma", positive=True)

    >>> X = LogNormal("x", mu, sigma)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /                         2\
          |          -(-mu + log(x)) |
          |          ----------------|
          |                     2    |
          |     ___      2*sigma     |
          |   \/ 2 *e                |
    Lambda|x, -----------------------|
          |             ____         |
          \       2*x*\/ pi *sigma   /

    >>> X = LogNormal('x', 0, 1) # Mean 0, standard deviation 1

    >>> density(X)
    Lambda(_x, sqrt(2)*exp(-log(_x)**2/2)/(2*_x*sqrt(pi)))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Lognormal
    [2] http://mathworld.wolfram.com/LogNormalDistribution.html
    """

    return LogNormalPSpace(name, mean, std).value

#-------------------------------------------------------------------------------
# Maxwell distribution ---------------------------------------------------------

class MaxwellPSpace(SingleContinuousPSpace):
    def __new__(cls, name, a):
        a = sympify(a)

        x = Symbol(name)

        pdf = sqrt(2/pi)*x**2*exp(-x**2/(2*a**2))/a**3
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Maxwell(name, a):
    r"""
    Create a continuous random variable with a Maxwell distribution.

    The density of the Maxwell distribution is given by

    .. math::
        f(x) := \sqrt{\frac{2}{\pi}} \frac{x^2 e^{-x^2/(2a^2)}}{a^3}

    with :math:`x \geq 0`.

    Parameters
    ==========

    a : Real number, `a` > 0

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Maxwell, density, E, variance
    >>> from sympy import Symbol, simplify

    >>> a = Symbol("a", positive=True)

    >>> X = Maxwell("x", a)

    >>> density(X)
    Lambda(_x, sqrt(2)*_x**2*exp(-_x**2/(2*a**2))/(sqrt(pi)*a**3))

    >>> E(X)
    2*sqrt(2)*a/sqrt(pi)

    >>> simplify(variance(X))
    a**2*(-8 + 3*pi)/pi

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Maxwell_distribution
    [2] http://mathworld.wolfram.com/MaxwellDistribution.html
    """

    return MaxwellPSpace(name, a).value

#-------------------------------------------------------------------------------
# Nakagami distribution --------------------------------------------------------

class NakagamiPSpace(SingleContinuousPSpace):
    def __new__(cls, name, mu, omega):
        mu, omega = sympify(mu), sympify(omega)
        x = Symbol(name)
        pdf = 2*mu**mu/(gamma(mu)*omega**mu)*x**(2*mu-1)*exp(-mu/omega*x**2)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Nakagami(name, mu, omega):
    r"""
    Create a continuous random variable with a Nakagami distribution.

    The density of the Nakagami distribution is given by

    .. math::
        f(x) := \frac{2\mu^\mu}{\Gamma(\mu)\omega^\mu} x^{2\mu-1}
                \exp\left(-\frac{\mu}{\omega}x^2 \right)

    with :math:`x > 0`.

    Parameters
    ==========

    mu : Real number, :math:`mu \geq \frac{1}{2}` a shape
    omega : Real number, `omega` > 0 the spread

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Nakagami, density, E, variance
    >>> from sympy import Symbol, simplify, pprint

    >>> mu = Symbol("mu", positive=True)
    >>> omega = Symbol("omega", positive=True)

    >>> X = Nakagami("x", mu, omega)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /                                2   \
          |                              -x *mu|
          |                              ------|
          |      2*mu - 1   mu      -mu  omega |
          |   2*x        *mu  *omega   *e      |
    Lambda|x, ---------------------------------|
          \               gamma(mu)            /

    >>> simplify(E(X, meijerg=True))
    sqrt(mu)*sqrt(omega)*gamma(mu + 1/2)/gamma(mu + 1)

    >>> V = simplify(variance(X, meijerg=True))
    >>> pprint(V, use_unicode=False)
          /                               2          \
    omega*\gamma(mu)*gamma(mu + 1) - gamma (mu + 1/2)/
    --------------------------------------------------
                 gamma(mu)*gamma(mu + 1)

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Nakagami_distribution
    """

    return NakagamiPSpace(name, mu, omega).value

#-------------------------------------------------------------------------------
# Normal distribution ----------------------------------------------------------

class NormalPSpace(SingleContinuousPSpace):
    def __new__(cls, name, mean, std):
        mean, std = sympify(mean), sympify(std)

        _value_check(std > 0, "Standard deviation must be positive")

        x = Symbol(name)
        pdf = exp(-(x-mean)**2 / (2*std**2)) / (sqrt(2*pi)*std)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.std = std
        obj.variance = std**2
        return obj

    def sample(self):
        return {self.value: random.normalvariate(self.mean, self.std)}

def Normal(name, mean, std):
    r"""
    Create a continuous random variable with a Normal distribution.

    The density of the Normal distribution is given by

    .. math::
        f(x) := \frac{1}{\sigma\sqrt{2\pi}} e^{ -\frac{(x-\mu)^2}{2\sigma^2} }

    Parameters
    ==========

    mu : Real number, the mean
    sigma : Real number, :math:`\sigma^2 > 0` the variance

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Normal, density, E, std, cdf, skewness
    >>> from sympy import Symbol, simplify, pprint

    >>> mu = Symbol("mu")
    >>> sigma = Symbol("sigma", positive=True)

    >>> X = Normal("x", mu, sigma)

    >>> density(X)
    Lambda(_x, sqrt(2)*exp(-(_x - mu)**2/(2*sigma**2))/(2*sqrt(pi)*sigma))

    >>> C = simplify(cdf(X))
    >>> pprint(C, use_unicode=False)
              /      /  ___         \    \
              |      |\/ 2 *(z - mu)|    |
              |   erf|--------------|    |
              |      \   2*sigma    /   1|
        Lambda|z, ------------------- + -|
              \            2            2/

    >>> simplify(skewness(X))
    0

    >>> X = Normal("x", 0, 1) # Mean 0, standard deviation 1
    >>> density(X)
    Lambda(_x, sqrt(2)*exp(-_x**2/2)/(2*sqrt(pi)))

    >>> E(2*X + 1)
    1

    >>> simplify(std(2*X + 1))
    2

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Normal_distribution
    [2] http://mathworld.wolfram.com/NormalDistributionFunction.html
    """

    return NormalPSpace(name, mean, std).value

#-------------------------------------------------------------------------------
# Pareto distribution ----------------------------------------------------------

class ParetoPSpace(SingleContinuousPSpace):
    def __new__(cls, name, xm, alpha):
        xm, alpha = sympify(xm), sympify(alpha)

        _value_check(xm > 0, "Xm must be positive")
        _value_check(alpha > 0, "Alpha must be positive")

        x = Symbol(name)
        pdf = alpha * xm**alpha / x**(alpha+1)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(xm, oo))
        obj.xm = xm
        obj.alpha = alpha
        return obj

    def sample(self):
        return {self.value: random.paretovariate(self.alpha)}

def Pareto(name, xm, alpha):
    r"""
    Create a continuous random variable with the Pareto distribution.

    The density of the Pareto distribution is given by

    .. math::
        f(x) := \frac{\alpha\,x_\mathrm{m}^\alpha}{x^{\alpha+1}}

    with :math:`x \in [x_m,\infty]`.

    Parameters
    ==========

    xm : Real number, `xm` > 0 a scale
    alpha : Real number, `alpha` > 0 a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Pareto, density
    >>> from sympy import Symbol

    >>> xm = Symbol("xm", positive=True)
    >>> beta = Symbol("beta", positive=True)

    >>> X = Pareto("x", xm, beta)

    >>> density(X)
    Lambda(_x, _x**(-beta - 1)*beta*xm**beta)

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Pareto_distribution
    [2] http://mathworld.wolfram.com/ParetoDistribution.html
    """

    return ParetoPSpace(name, xm, alpha).value

#-------------------------------------------------------------------------------
# Rayleigh distribution --------------------------------------------------------

class RayleighPSpace(SingleContinuousPSpace):
    def __new__(cls, name, sigma):
        sigma = sympify(sigma)
        x = Symbol(name)
        pdf = x/sigma**2*exp(-x**2/(2*sigma**2))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Rayleigh(name, sigma):
    r"""
    Create a continuous random variable with a Rayleigh distribution.

    The density of the Rayleigh distribution is given by

    .. math ::
        f(x) := \frac{x}{\sigma^2} e^{-x^2/2\sigma^2}

    with :math:`x > 0`.

    Parameters
    ==========

    sigma : Real number, `sigma` > 0

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Rayleigh, density, E, variance
    >>> from sympy import Symbol, simplify

    >>> sigma = Symbol("sigma", positive=True)

    >>> X = Rayleigh("x", sigma)

    >>> density(X)
    Lambda(_x, _x*exp(-_x**2/(2*sigma**2))/sigma**2)

    >>> E(X)
    sqrt(2)*sqrt(pi)*sigma/2

    >>> variance(X)
    -pi*sigma**2/2 + 2*sigma**2

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Rayleigh_distribution
    [2] http://mathworld.wolfram.com/RayleighDistribution.html
    """

    return RayleighPSpace(name, sigma).value

#-------------------------------------------------------------------------------
# StudentT distribution --------------------------------------------------------

class StudentTPSpace(SingleContinuousPSpace):
    def __new__(cls, name, nu):
        nu = sympify(nu)
        x = Symbol(name)
        pdf = 1/(sqrt(nu)*beta_fn(S(1)/2,nu/2))*(1+x**2/nu)**(-(nu+1)/2)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def StudentT(name, nu):
    r"""
    Create a continuous random variable with a student's t distribution.

    The density of the student's t distribution is given by

    .. math::
        f(x) := \frac{\Gamma \left(\frac{\nu+1}{2} \right)}
                {\sqrt{\nu\pi}\Gamma \left(\frac{\nu}{2} \right)}
                \left(1+\frac{x^2}{\nu} \right)^{-\frac{\nu+1}{2}}

    Parameters
    ==========

    nu : Real number, `nu` > 0, the degrees of freedom

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import StudentT, density, E, variance
    >>> from sympy import Symbol, simplify, pprint

    >>> nu = Symbol("nu", positive=True)

    >>> X = StudentT("x", nu)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /             nu   1              \
          |           - -- - -              |
          |             2    2              |
          |   / 2    \                      |
          |   |x     |              /nu   1\|
          |   |-- + 1|        *gamma|-- + -||
          |   \nu    /              \2    2/|
    Lambda|x, ------------------------------|
          |        ____   ____      /nu\    |
          |      \/ pi *\/ nu *gamma|--|    |
          \                         \2 /    /

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Student_t-distribution
    [2] http://mathworld.wolfram.com/Studentst-Distribution.html
    """

    return StudentTPSpace(name, nu).value

#-------------------------------------------------------------------------------
# Triangular distribution ------------------------------------------------------

class TriangularPSpace(SingleContinuousPSpace):
    def __new__(cls, name, a, b, c):
        a, b, c = sympify(a), sympify(b), sympify(c)

        x = Symbol(name)
        pdf = Piecewise(
                (2*(x-a)/((b-a)*(c-a)), And(a<=x, x<c)),
                (2/(b-a), Eq(x,c)),
                (2*(b-x)/((b-a)*(b-c)), And(c<x, x<=b)),
                (S.Zero, True))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Triangular(name, a, b, c):
    r"""
    Create a continuous random variable with a triangular distribution.

    The density of the triangular distribution is given by

    .. math::
        f(x) := \begin{cases}
                  0 & \mathrm{for\ } x < a, \\
                  \frac{2(x-a)}{(b-a)(c-a)} & \mathrm{for\ } a \le x < c, \\
                  \frac{2}{b-a} & \mathrm{for\ } x = c, \\
                  \frac{2(b-x)}{(b-a)(b-c)} & \mathrm{for\ } c < x \le b, \\
                  0 & \mathrm{for\ } b < x.
                \end{cases}

    Parameters
    ==========

    a : Real number, :math:`a \in \left(-\infty, \infty\right)`
    b : Real number, :math:`a < b`
    c : Real number, :math:`a \leq c \leq b`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Triangular, density, E
    >>> from sympy import Symbol

    >>> a = Symbol("a")
    >>> b = Symbol("b")
    >>> c = Symbol("c")

    >>> X = Triangular("x", a,b,c)

    >>> density(X)
    Lambda(_x, Piecewise(((2*_x - 2*a)/((-a + b)*(-a + c)),
                        And(_x < c, a <= _x)),
                        (2/(-a + b), _x == c),
                        ((-2*_x + 2*b)/((-a + b)*(b - c)),
                        And(_x <= b, c < _x)),
                        (0, True)))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Triangular_distribution
    [2] http://mathworld.wolfram.com/TriangularDistribution.html
    """

    return TriangularPSpace(name, a, b, c).value

#-------------------------------------------------------------------------------
# Uniform distribution ---------------------------------------------------------

class UniformPSpace(SingleContinuousPSpace):
    def __new__(cls, name, left, right):
        left, right = sympify(left), sympify(right)

        x = Symbol(name)
        pdf = Piecewise(
                (S.One/(right-left), And(left <= x, x <= right)),
                (S.Zero, True))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.left = left
        obj.right = right
        return obj

    def compute_cdf(self, expr, **kwargs):
        from sympy import Lambda, Min
        z = Dummy('z', real=True, bounded=True)
        result = SingleContinuousPSpace.compute_cdf(self, expr, **kwargs)
        result = result(z).subs({Min(z, self.right): z,
                                 Min(z, self.left, self.right): self.left})
        return Lambda(z, result)

    def integrate(self, expr, rvs=None, **kwargs):
        from sympy import Max, Min
        result = SingleContinuousPSpace.integrate(self, expr, rvs, **kwargs)
        result = result.subs({Max(self.left, self.right): self.right,
                              Min(self.left, self.right): self.left})
        return result

    def sample(self):
        return {self.value: random.uniform(self.left, self.right)}

def Uniform(name, left, right):
    r"""
    Create a continuous random variable with a uniform distribution.

    The density of the uniform distribution is given by

    .. math::
        f(x) := \begin{cases}
                  \frac{1}{b - a} & \text{for } x \in [a,b]  \\
                  0               & \text{otherwise}
                \end{cases}

    with :math:`x \in [a,b]`.

    Parameters
    ==========

    a : Real number, :math:`-\infty < a` the left boundary
    b : Real number, :math:`a < b < \infty` the right boundary

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Uniform, density, cdf, E, variance, skewness
    >>> from sympy import Symbol, simplify

    >>> a = Symbol("a")
    >>> b = Symbol("b")

    >>> X = Uniform("x", a, b)

    >>> density(X)
    Lambda(_x, Piecewise((1/(-a + b), And(_x <= b, a <= _x)), (0, True)))

    >>> cdf(X)
    Lambda(_z, _z/(-a + b) - a/(-a + b))

    >>> simplify(E(X))
    a/2 + b/2

    >>> simplify(variance(X))
    a**2/12 - a*b/6 + b**2/12

    >>> simplify(skewness(X))
    0

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Uniform_distribution_%28continuous%29
    [2] http://mathworld.wolfram.com/UniformDistribution.html
    """

    return UniformPSpace(name, left, right).value

#-------------------------------------------------------------------------------
# UniformSum distribution ------------------------------------------------------

class UniformSumPSpace(SingleContinuousPSpace):
    def __new__(cls, name, n):
        n = sympify(n)

        x = Symbol(name)
        k = Dummy("k")
        pdf =1/factorial(n-1)*Sum((-1)**k*binomial(n,k)*(x-k)**(n-1), (k,0,floor(x)))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0,n))
        return obj

def UniformSum(name, n):
    r"""
    Create a continuous random variable with an Irwin-Hall distribution.

    The probability distribution function depends on a single parameter
    `n` which is an integer.

    The density of the Irwin-Hall distribution is given by

    .. math ::
        f(x) := \frac{1}{(n-1)!}\sum_{k=0}^{\lfloor x\rfloor}(-1)^k
                \binom{n}{k}(x-k)^{n-1}

    Parameters
    ==========

    n : Integral number, `n` > 0

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import UniformSum, density
    >>> from sympy import Symbol, pprint

    >>> n = Symbol("n", integer=True)

    >>> X = UniformSum("x", n)

    >>> D = density(X)
    >>> pprint(D, use_unicode=False)
          /   floor(x)                        \
          |     ___                           |
          |     \  `                          |
          |      \         k         n - 1 /n\|
          |       )    (-1) *(-k + x)     *| ||
          |      /                         \k/|
          |     /__,                          |
          |    k = 0                          |
    Lambda|x, --------------------------------|
          \               (n - 1)!            /

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Uniform_sum_distribution
    [2] http://mathworld.wolfram.com/UniformSumDistribution.html
    """

    return UniformSumPSpace(name, n).value

#-------------------------------------------------------------------------------
# Weibull distribution ---------------------------------------------------------

class WeibullPSpace(SingleContinuousPSpace):
    def __new__(cls, name, alpha, beta):
        alpha, beta = sympify(alpha), sympify(beta)

        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(beta > 0, "Beta must be positive")

        x = Symbol(name)
        pdf = beta * (x/alpha)**(beta-1) * exp(-(x/alpha)**beta) / alpha

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.weibullvariate(self.alpha, self.beta)}

def Weibull(name, alpha, beta):
    r"""
    Create a continuous random variable with a Weibull distribution.

    The density of the Weibull distribution is given by

    .. math::
        f(x) := \begin{cases}
                  \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}
                  e^{-(x/\lambda)^{k}} & x\geq0\\
                  0 & x<0
                \end{cases}

    Parameters
    ==========

    lambda : Real number, :math:`\lambda > 0` a scale
    k : Real number, `k` > 0 a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Weibull, density, E, variance
    >>> from sympy import Symbol, simplify

    >>> l = Symbol("lambda", positive=True)
    >>> k = Symbol("k", positive=True)

    >>> X = Weibull("x", l, k)

    >>> density(X)
    Lambda(_x, k*(_x/lambda)**(k - 1)*exp(-(_x/lambda)**k)/lambda)

    >>> simplify(E(X))
    lambda*gamma(1 + 1/k)

    >>> simplify(variance(X))
    lambda**2*(-gamma(1 + 1/k)**2 + gamma(1 + 2/k))

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Weibull_distribution
    [2] http://mathworld.wolfram.com/WeibullDistribution.html

    """

    return WeibullPSpace(name, alpha, beta).value

#-------------------------------------------------------------------------------
# Wigner semicircle distribution -----------------------------------------------

class WignerSemicirclePSpace(SingleContinuousPSpace):
    def __new__(cls, name, R):
        R = sympify(R)

        x = Symbol(name)
        pdf = 2/(pi*R**2)*sqrt(R**2-x**2)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(-R, R))
        return obj

def WignerSemicircle(name, R):
    r"""
    Create a continuous random variable with a Wigner semicircle distribution.

    The density of the Wigner semicircle distribution is given by

    .. math::
        f(x) := \frac2{\pi R^2}\,\sqrt{R^2-x^2}

    with :math:`x \in [-R,R]`.

    Parameters
    ==========

    R : Real number, `R` > 0 the radius

    Returns
    =======

    A `RandomSymbol`.

    Examples
    ========

    >>> from sympy.stats import WignerSemicircle, density, E
    >>> from sympy import Symbol, simplify

    >>> R = Symbol("R", positive=True)

    >>> X = WignerSemicircle("x", R)

    >>> density(X)
    Lambda(_x, 2*sqrt(-_x**2 + R**2)/(pi*R**2))

    >>> E(X)
    0

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    [2] http://mathworld.wolfram.com/WignersSemicircleLaw.html
    """

    return WignerSemicirclePSpace(name, R).value
