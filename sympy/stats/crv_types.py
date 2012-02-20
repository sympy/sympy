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
                   Piecewise, And, Eq, binomial, factorial, Sum, floor, Abs, log)
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

    >>> x = Symbol('x')
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
    def __new__(cls, a, b, symbol=None):
        a, b = sympify(a), sympify(b)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 1/(pi*sqrt((x-a)*(b-x)))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(a, b))
        return obj

def Arcsin(a=0, b=1, symbol=None):
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

    A `RandomSymbol` X.

    Examples
    ========

    >>> from sympy.stats import Arcsin, Density
    >>> from sympy import Symbol, simplify

    >>> a = Symbol("a", real=True)
    >>> b = Symbol("b", real=True)
    >>> x = Symbol("x")

    >>> X = Arcsin(a, b, symbol=x)

    >>> Density(X)
    (x, 1/(pi*sqrt((-a + x)*(b - x))))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Arcsine_distribution
    """

    return ArcsinPSpace(a, b, symbol).value

#-------------------------------------------------------------------------------
# Benini distribution ----------------------------------------------------------

class BeniniPSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, sigma, symbol = None):
        alpha, beta, sigma = sympify(alpha), sympify(beta), sympify(sigma)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-alpha*log(x/sigma)-beta*log(x/sigma)**2)*(alpha/x+2*beta*log(x/sigma)/x)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(sigma, oo))
        return obj

def Benini(alpha, beta, sigma, symbol=None):
    r"""
    Create a Continuous Random Variable with a Benini distribution.

    The density of the Benini distribution is given by

    .. math::
        f(x) := e^{-\alpha\log{\frac{x}{\sigma}}-\beta\log\left[{\frac{x}{\sigma}}\right]^2}
        \left(\frac{\alpha}{x}+\frac{2\beta\log{\frac{x}{\sigma}}}{x}\right)

    Parameters
    ==========

    alpha : Real number, `alpha` > 0 a shape
    beta : Real number, `beta` > 0 a shape
    sigma : Real number, `sigma` > 0 a scale

    Returns
    =======

    A `RandomSymbol` X.

    Examples
    ========

    >>> from sympy.stats import Benini, Density
    >>> from sympy import Symbol, simplify

    >>> alpha = Symbol("alpha", positive=True)
    >>> beta = Symbol("beta", positive=True)
    >>> sigma = Symbol("sigma", positive=True)
    >>> x = Symbol("x")

    >>> X = Benini(alpha, beta, sigma, symbol=x)

    >>> Density(X)
    (x, (alpha/x + 2*beta*log(x/sigma)/x)*exp(-alpha*log(x/sigma) - beta*log(x/sigma)**2))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Benini_distribution
    """

    return BeniniPSpace(alpha, beta, sigma, symbol).value

#-------------------------------------------------------------------------------
# Beta distribution ------------------------------------------------------------

class BetaPSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, symbol=None):
        alpha, beta = sympify(alpha), sympify(beta)

        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(beta > 0, "Beta must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(alpha-1) * (1-x)**(beta-1) / beta_fn(alpha, beta)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, 1))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.betavariate(self.alpha, self.beta)}

def Beta(alpha, beta, symbol=None):
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

    A `RandomSymbol` X.

    Examples
    ========

    >>> from sympy.stats import Beta, Density, E, Var
    >>> from sympy import Symbol, simplify

    >>> alpha = Symbol("alpha", positive=True)
    >>> beta = Symbol("beta", positive=True)
    >>> x = Symbol("x")

    >>> X = Beta(alpha, beta, symbol=x)

    >>> Density(X)
    Lambda(_x, _x**(alpha - 1)*(-_x + 1)**(beta - 1)*gamma(alpha + beta)/(gamma(alpha)*gamma(beta)))

    >>> simplify(E(X, meijerg=True))
    alpha/(alpha + beta)

    >>> simplify(Var(X, meijerg=True))
    alpha*beta/((alpha + beta)**2*(alpha + beta + 1))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Beta_distribution
    .. [2] http://mathworld.wolfram.com/BetaDistribution.html
    """

    return BetaPSpace(alpha, beta, symbol).value

#-------------------------------------------------------------------------------
# Beta prime distribution ------------------------------------------------------

class BetaPrimePSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, symbol=None):
        alpha, beta = sympify(alpha), sympify(beta)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(alpha-1)*(1+x)**(-alpha-beta)/beta_fn(alpha, beta)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def BetaPrime(alpha, beta, symbol=None):
    r"""
    Create a Continuous Random Variable with a Beta prime distribution.

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

    A `RandomSymbol` X.

    Examples
    ========

    >>> from sympy.stats import BetaPrime, Density
    >>> from sympy import Symbol

    >>> alpha = Symbol("alpha", positive=True)
    >>> beta = Symbol("beta", positive=True)
    >>> x = Symbol("x")

    >>> X = BetaPrime(alpha, beta, symbol=x)

    >>> Density(X)
    (x, x**(alpha - 1)*(x + 1)**(-alpha - beta)*gamma(alpha + beta)/(gamma(alpha)*gamma(beta)))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Beta_prime_distribution
    .. [2] http://mathworld.wolfram.com/BetaPrimeDistribution.html
    """

    return BetaPrimePSpace(alpha, beta, symbol).value

#-------------------------------------------------------------------------------
# Cauchy distribution ----------------------------------------------------------

class CauchyPSpace(SingleContinuousPSpace):
    def __new__(cls, x0, gamma, symbol = None):
        x0, gamma = sympify(x0), sympify(gamma)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 1/(pi*gamma*(1+((x-x0)/gamma)**2))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Cauchy(x0, gamma, symbol=None):
    r"""
    Create a Continuous Random Variable with a Cauchy distribution.

    The density of the Cauchy distribution is given by

    .. math::
        f(x) := \frac{1}{\pi} \arctan\left(\frac{x-x_0}{\gamma}\right)+\frac{1}{2}

    Parameters
    ==========

    x0 : Real number, the location
    gamma : Real number, `gamma` > 0 the scale

    Returns
    =======

    A `RandomSymbol` X.

    Examples
    ========

    >>> from sympy.stats import Cauchy, Density
    >>> from sympy import Symbol

    >>> x0 = Symbol("x0")
    >>> gamma = Symbol("gamma", positive=True)
    >>> x = Symbol("x")

    >>> X = Cauchy(x0, gamma, symbol=x)

    >>> Density(X)
    (x, 1/(pi*gamma*(1 + (x - x0)**2/gamma**2)))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Cauchy_distribution
    .. [2] http://mathworld.wolfram.com/CauchyDistribution.html
    """

    return CauchyPSpace(x0, gamma, symbol).value

#-------------------------------------------------------------------------------
# Chi distribution -------------------------------------------------------------

class ChiPSpace(SingleContinuousPSpace):
    def __new__(cls, k, symbol = None):
        k = sympify(k)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 2**(1-k/2)*x**(k-1)*exp(-x**2/2)/gamma(k/2)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Chi(k, symbol=None):
    r"""
    Create a Continuous Random Variable with a Chi distribution.

    The density of the Chi distribution is given by

    .. math::
        f(x) := \frac{2^{1-k/2}x^{k-1}e^{-x^2/2}}{\Gamma(k/2)}

    with :math:`x \geq 0`.

    Parameters
    ==========

    k : Integer, `k` > 0 the number of degrees of freedom

    Returns
    =======

    A `RandomSymbol` X.

    Examples
    ========

    >>> from sympy.stats import Chi, Density, E, Std
    >>> from sympy import Symbol, simplify

    >>> k = Symbol("k", integer=True)
    >>> x = Symbol("x")

    >>> X = Chi(k, symbol=x)

    >>> Density(X)
    (x, 2**(-k/2 + 1)*x**(k - 1)*exp(-x**2/2)/gamma(k/2))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Chi_distribution
    .. [2] http://mathworld.wolfram.com/ChiDistribution.html
    """

    return ChiPSpace(k, symbol).value

#-------------------------------------------------------------------------------
# Dagum distribution -----------------------------------------------------

class DagumPSpace(SingleContinuousPSpace):
    def __new__(cls, p, a, b, symbol = None):
        p, a, b = sympify(p), sympify(a), sympify(b)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = a*p/x*((x/b)**(a*p)/(((x/b)**a+1)**(p+1)))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Dagum(p, a, b, symbol=None):
    """
    Create a Continuous Random Variable with a Dagum distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Dagum, Density, E, Std
    >>> from sympy import Symbol, simplify
    """

    return DagumPSpace(p, a, b, symbol).value

#-------------------------------------------------------------------------------
# Exponential distribution -----------------------------------------------------

class ExponentialPSpace(SingleContinuousPSpace):
    def __new__(cls, rate, symbol=None):
        rate = sympify(rate)

        _value_check(rate > 0, "Rate must be positive.")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = rate * exp(-rate*x)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.rate = rate
        return obj

    def sample(self):
        return {self.value: random.expovariate(self.rate)}

def Exponential(rate, symbol=None):
    """
    Create a Continuous Random Variable with an Exponential distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Exponential, Density, E, Std
    >>> from sympy import Symbol

    >>> X = Exponential(rate=10, symbol=Symbol('x')) # Decay rate equals 10
    >>> Density(X)
    Lambda(_x, 10*exp(-10*_x))

    >>> E(X)
    1/10

    >>> Std(X)
    1/10
    """

    return ExponentialPSpace(rate, symbol).value

#-------------------------------------------------------------------------------
# Gamma distribution -----------------------------------------------------------

class GammaPSpace(SingleContinuousPSpace):
    def __new__(cls, k, theta, symbol=None):
        k, theta = sympify(k), sympify(theta)

        _value_check(k > 0, "k must be positive")
        _value_check(theta > 0, "Theta must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(k-1) * exp(-x/theta) / (gamma(k)*theta**k)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.k = k
        obj.theta = theta
        return obj

    def sample(self):
        return {self.value: random.gammavariate(self.k, self.theta)}

def Gamma(k, theta, symbol=None):
    """
    Create a Continuous Random Variable with a Gamma distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Gamma, Density, E, Std
    >>> from sympy import symbols
    >>> x, k, theta = symbols('x k theta', positive=True)

    >>> X = Gamma(k, theta, symbol=x)
    >>> Density(X)
    Lambda(_x, _x**(k - 1)*theta**(-k)*exp(-_x/theta)/gamma(k))

    >>> E(X)
    theta*gamma(k + 1)/gamma(k)
    """

    return GammaPSpace(k, theta, symbol).value

#-------------------------------------------------------------------------------
# Laplace distribution ---------------------------------------------------------

class LaplacePSpace(SingleContinuousPSpace):
    def __new__(cls, mu, b, symbol=None):
        mu, b = sympify(mu), sympify(b)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 1/(2*b)*exp(-Abs(x-mu)/b)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Laplace(mu, b, symbol=None):
    """
    Create a Continuous Random Variable with a Laplace distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Laplace, Density, E, Std
    >>> from sympy import Symbol
    """

    return LaplacePSpace(mu, b, symbol).value

#-------------------------------------------------------------------------------
# Logistic distribution --------------------------------------------------------

class LogisticPSpace(SingleContinuousPSpace):
    def __new__(cls, mu, s, symbol=None):
        mu, s = sympify(mu), sympify(s)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-(x-mu)/s)/(s*(1+exp(-(x-mu)/s))**2)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Logistic(mu, s, symbol=None):
    """
    Create a Continuous Random Variable with a Logistic distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Logistic, Density, E, Std
    >>> from sympy import Symbol
    """

    return LogisticPSpace(mu, s, symbol).value

#-------------------------------------------------------------------------------
# Log Normal distribution ------------------------------------------------------

class LogNormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, std, symbol=None):
        mean, std = sympify(mean), sympify(std)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-(log(x)-mean)**2 / (2*std**2)) / (x*sqrt(2*pi)*std)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.mean = mean
        obj.std = std
        return obj

    def sample(self):
        return {self.value: random.lognormvariate(self.mean, self.std)}

def LogNormal(mean, std, symbol=None):
    """
    Create a Continuous Random Variable with a LogNormal distribution.

    Note: Only density and sampling work.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import LogNormal, Density, E, Std
    >>> from sympy import Symbol, simplify

    >>> X = LogNormal(0, 1, symbol=Symbol('x')) # Mean 0, standard deviation 1

    >>> Density(X)
    Lambda(_x, sqrt(2)*exp(-log(_x)**2/2)/(2*_x*sqrt(pi)))
    """

    return LogNormalPSpace(mean, std, symbol).value

#-------------------------------------------------------------------------------
# Maxwell distribution ---------------------------------------------------------

class MaxwellPSpace(SingleContinuousPSpace):
    def __new__(cls, a, symbol = None):
        a = sympify(a)

        x = symbol or SingleContinuousPSpace.create_symbol()

        pdf = sqrt(2/pi)*x**2*exp(-x**2/(2*a**2))/a**3
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Maxwell(a, symbol=None):
    """
    Create a Continuous Random Variable with a Maxwell distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Maxwell, Density, E, Std
    >>> from sympy import Symbol, simplify
    """

    return MaxwellPSpace(a, symbol).value

#-------------------------------------------------------------------------------
# Nakagami distribution --------------------------------------------------------

class NakagamiPSpace(SingleContinuousPSpace):
    def __new__(cls, mu, omega, symbol = None):
        mu, omega = sympify(mu), sympify(omega)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 2*mu**mu/(gamma(mu)*omega**mu)*x**(2*mu-1)*exp(-mu/omega*x**2)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Nakagami(mu, omega, symbol=None):
    """
    Create a Continuous Random Variable with a Nakagami distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Nakagami, Density, E, Std
    >>> from sympy import Symbol, simplify
    """

    return NakagamiPSpace(mu, omega, symbol).value

#-------------------------------------------------------------------------------
# Normal distribution ----------------------------------------------------------

class NormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, std, symbol=None):
        mean, std = sympify(mean), sympify(std)

        _value_check(std > 0, "Standard deviation must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-(x-mean)**2 / (2*std**2)) / (sqrt(2*pi)*std)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.std = std
        obj.variance = std**2
        return obj

    def sample(self):
        return {self.value: random.normalvariate(self.mean, self.std)}

def Normal(mean, std, symbol=None):
    """
    Create a Continuous Random Variable with a Normal distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Normal, Density, E, Std
    >>> from sympy import Symbol, simplify

    >>> X = Normal(0, 1, symbol=Symbol('x')) # Mean 0, standard deviation 1
    >>> Density(X)
    Lambda(_x, sqrt(2)*exp(-_x**2/2)/(2*sqrt(pi)))

    >>> E(2*X + 1)
    1

    >>> simplify(Std(2*X + 1))
    2
    """

    return NormalPSpace(mean, std, symbol).value

#-------------------------------------------------------------------------------
# Pareto distribution ----------------------------------------------------------

class ParetoPSpace(SingleContinuousPSpace):
    def __new__(cls, xm, alpha, symbol=None):
        xm, alpha = sympify(xm), sympify(alpha)

        _value_check(xm > 0, "Xm must be positive")
        _value_check(alpha > 0, "Alpha must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = alpha * xm**alpha / x**(alpha+1)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(xm, oo))
        obj.xm = xm
        obj.alpha = alpha
        return obj

    def sample(self):
        return {self.value: random.paretovariate(self.alpha)}

def Pareto(xm, alpha, symbol=None):
    """
    Create a Continuous Random Variable with the Pareto distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Pareto, Density, E, Std
    >>> from sympy import symbols

    >>> x, xm, beta = symbols('x xm beta', positive=True)
    >>> X = Pareto(xm, beta, symbol=x)
    >>> Density(X)
    Lambda(_x, _x**(-beta - 1)*beta*xm**beta)
    """

    return ParetoPSpace(xm, alpha, symbol).value

#-------------------------------------------------------------------------------
# Rayleigh distribution --------------------------------------------------------

class RayleighPSpace(SingleContinuousPSpace):
    def __new__(cls, sigma, symbol = None):
        sigma = sympify(sigma)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x/sigma**2*exp(-x**2/(2*sigma**2))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        return obj

def Rayleigh(sigma, symbol=None):
    """
    Create a Continuous Random Variable with a Rayleigh distribution.

    The density of the Rayleigh distribution is given by

    .. math ::
        \frac{x}{\sigma^2} e^{-x^2/2\sigma^2}.

    Parameters
    ==========

    sigma : Real number, `sigma` > 0

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Rayleigh, Density, E, Var
    >>> from sympy import Symbol, simplify

    >>> sigma = Symbol("sigma", positive=True)
    >>> x = Symbol("x")

    >>> X = Rayleigh(sigma, symbol=x)

    >>> Density(X)
    (x, x*exp(-x**2/(2*sigma**2))/sigma**2)

    >>> E(X)
    sqrt(2)*sqrt(pi)*sigma/2

    >>> Var(X)
    -pi*sigma**2/2 + 2*sigma**2

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Rayleigh_distribution
    .. [2] http://mathworld.wolfram.com/RayleighDistribution.html
    """

    return RayleighPSpace(sigma, symbol).value

#-------------------------------------------------------------------------------
# StudentT distribution --------------------------------------------------------

class StudentTPSpace(SingleContinuousPSpace):
    def __new__(cls, nu, symbol = None):
        nu = sympify(nu)
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 1/(sqrt(nu)*beta_fn(S(1)/2,nu/2))*(1+x**2/nu)**(-(nu+1)/2)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def StudentT(nu, symbol=None):
    """
    Create a Continuous Random Variable with a student's t distribution.

    The density of the student's t distribution is given by

    .. math ::
        \frac{\Gamma \left(\frac{\nu+1}{2} \right)} {\sqrt{\nu\pi}\,\Gamma \left(\frac{\nu}{2} \right)}
        \left(1+\frac{x^2}{\nu} \right)^{-\frac{\nu+1}{2}}

    Parameters
    ==========

    nu : Real number, `nu` > 0
        Degrees of freedom

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import StudentT, Density, E, Var
    >>> from sympy import Symbol, simplify

    >>> nu = Symbol("nu", positive=True)
    >>> x = Symbol("x")

    >>> X = StudentT(nu, symbol=x)

    >>> Density(X)
    (x, (1 + x**2/nu)**(-nu/2 - 1/2)*gamma(nu/2 + 1/2)/(sqrt(pi)*sqrt(nu)*gamma(nu/2)))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Student_t
    .. [2] http://mathworld.wolfram.com/Studentst-Distribution.html
    """

    return StudentTPSpace(nu, symbol).value

#-------------------------------------------------------------------------------
# Triangular distribution ------------------------------------------------------

class TriangularPSpace(SingleContinuousPSpace):
    def __new__(cls, a, b, c, symbol=None):
        a, b, c = sympify(a), sympify(b), sympify(c)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = Piecewise(
                (2*(x-a)/((b-a)*(c-a)), And(a<=x, x<c)),
                (2/(b-a), Eq(x,c)),
                (2*(b-x)/((b-a)*(b-c)), And(c<x, x<=b)),
                (S.Zero, True))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        return obj

def Triangular(a, b, c, symbol=None):
    """
    Create a Continuous Random Variable with a Triangular distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Triangular, Density, E, Var
    >>> from sympy import symbols, simplify
    >>> x, a, b, c = symbols('x a b c')

    >>> X = Triangular(a, b, c, symbol=x)
    """

    return TriangularPSpace(a, b, c, symbol).value

#-------------------------------------------------------------------------------
# Uniform distribution ---------------------------------------------------------

class UniformPSpace(SingleContinuousPSpace):
    def __new__(cls, left, right, symbol=None):
        left, right = sympify(left), sympify(right)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = Piecewise(
                (S.Zero, x<left),
                (S.Zero, x>right),
                (S.One/(right-left), True))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.left = left
        obj.right = right
        return obj

    def sample(self):
        return {self.value: random.uniform(self.left, self.right)}

def Uniform(left, right, symbol=None):
    """
    Create a Continuous Random Variable with a Uniform distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Uniform, Density, E, Var
    >>> from sympy import symbols, simplify
    >>> x, l, r = symbols('x l r')

    >>> X = Uniform(l, r, symbol=x)

    >>> Density(X)
    Lambda(_x, Piecewise((0, _x < l), (0, _x > r), (1/(-l + r), True)))

    >>> simplify(E(X))
    l/2 + r/2

    >>> simplify(Var(X))
    l**2/12 - l*r/6 + r**2/12
    """

    return UniformPSpace(left, right, symbol).value

#-------------------------------------------------------------------------------
# UniformSum distribution ------------------------------------------------------

class UniformSumPSpace(SingleContinuousPSpace):
    def __new__(cls, n, symbol=None):
        n = sympify(n)

        x = symbol or SingleContinuousPSpace.create_symbol()
        k = Dummy("k")
        pdf =1/factorial(n-1)*Sum((-1)**k*binomial(n,k)*(x-k)**(n-1), (k,0,floor(x)))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0,n))
        return obj

def UniformSum(n, symbol=None):
    """
    Create a Continuous Random Variable with an Irwin-Hall distribution.

    The probability distribution function depends on a single parameter
    `n` which is an integer.

    The density is given by

    .. math ::
        \frac{1}{(n-1)!}\sum_{k=0}^{\lfloor x\rfloor}(-1)^k\binom{n}{k}(x-k)^{n-1}

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import UniformSum, Density, E, Var
    >>> from sympy import Symbol

    >>> n = Symbol("n", integer=True)
    >>> x = Symbol("x")

    >>> X = UniformSum(n, symbol=x)

    >>> Density(X)
    (x, Sum((-1)**_k*(-_k + x)**(n - 1)*binomial(n, _k), (_k, 0, floor(x)))/(n - 1)!)

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Uniform_sum_distribution
    """

    return UniformSumPSpace(n, symbol).value

#-------------------------------------------------------------------------------
# Weibull distribution ---------------------------------------------------------

class WeibullPSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, symbol=None):
        alpha, beta = sympify(alpha), sympify(beta)

        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(beta > 0, "Beta must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = beta * (x/alpha)**(beta-1) * exp(-(x/alpha)**beta) / alpha

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.weibullvariate(self.alpha, self.beta)}

def Weibull(alpha, beta, symbol=None):
    """
    Create a Continuous Random Variable with a Weibull distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Weibull, Density, E, Var
    >>> from sympy import symbols, simplify
    >>> x, a, b = symbols('x a b', positive=True)

    >>> X = Weibull(a, b, symbol=x)

    >>> Density(X)
    Lambda(_x, b*(_x/a)**(b - 1)*exp(-(_x/a)**b)/a)

    >>> simplify(E(X))
    a*gamma(1 + 1/b)

    >>> simplify(Var(X))
    -a**2*(gamma(1 + 1/b)**2 - gamma(1 + 2/b))
    """
    return WeibullPSpace(alpha, beta, symbol).value

#-------------------------------------------------------------------------------
# Wigner semicircle distribution -----------------------------------------------

class WignerSemicirclePSpace(SingleContinuousPSpace):
    def __new__(cls, R, symbol=None):
        R = sympify(R)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = 2/(pi*R**2)*sqrt(R**2-x**2)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(-R, R))
        return obj

def WignerSemicircle(R, symbol=None):
    """
    Create a Continuous Random Variable with a Wigner semicircle distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import WignerSemicircle, Density, E, Std
    >>> from sympy import Symbol, simplify
    """

    return WignerSemicirclePSpace(R, symbol).value
