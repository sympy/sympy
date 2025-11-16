"""
Alpha-stable distribution implementation for SymPy.
"""

from sympy.core.sympify import sympify
from sympy.stats.crv import SingleContinuousDistribution
from sympy.stats.crv_types import rv
from sympy.stats.rv import _value_check
from sympy import (exp, pi, tan, Abs, I, Piecewise, sqrt, oo,
                    Rational, Interval, log)


class AlphaStableDistribution(SingleContinuousDistribution):
    r"""
    Represents an alpha-stable distribution.

    Parameters
    ==========

    alpha : Real number, 0 < alpha <= 2
        Stability parameter (controls tail heaviness).
        When alpha < 2, the distribution has heavy tails.
    beta : Real number, -1 <= beta <= 1
        Skewness parameter. beta = 0 gives a symmetric distribution,
        beta > 0 skews right, beta < 0 skews left.
    scale : Positive real number
        Scale parameter (analogous to standard deviation for normal distribution).
    location : Real number
        Location parameter (analogous to mean for normal distribution).

    Notes
    =====

    The alpha-stable distribution is a family of continuous probability
    distributions that are closed under addition. They are characterized by
    four parameters: alpha (stability), beta (skewness), scale, and location.

    The distribution has no closed-form PDF for general parameters, but
    special cases exist:

    * alpha = 2, beta = 0: Normal distribution N(location, 2*scale^2)
    * alpha = 1, beta = 0: Cauchy distribution
    * alpha = 1/2, beta = 1: Levy distribution

    The characteristic function is available for all parameter values.

    Examples
    ========

    >>> from sympy.stats import AlphaStable, density, cdf
    >>> from sympy import Symbol, pprint
    >>> x = Symbol('x', real=True)

    Create a Cauchy distribution:

    >>> X = AlphaStable('X', 1, 0, 1, 0)
    >>> density(X)(x)
    1/(pi*(x**2 + 1))

    Create a standard normal distribution:

    >>> Y = AlphaStable('Y', 2, 0, 1, 0)
    >>> density(Y)(x)
    sqrt(2)*exp(-x**2/2)/(2*sqrt(pi))

    Create a general stable distribution:

    >>> Z = AlphaStable('Z', alpha=1.5, beta=0.5, scale=2, location=1)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Stable_distribution
    .. [2] Nolan, J. P. (2020). Univariate stable distributions.
           Springer Series in Operations Research and Financial Engineering.
    .. [3] Zolotarev, V. M. (1986). One-dimensional stable distributions.
           American Mathematical Society.

    """

    _argnames = ('alpha', 'beta', 'scale', 'location')

    @staticmethod
    def check(alpha, beta, scale, location):
        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(alpha <= 2, "Alpha must be less than or equal to 2")
        _value_check(beta >= -1, "Beta must be >= -1")
        _value_check(beta <= 1, "Beta must be <= 1")
        _value_check(scale > 0, "Scale must be positive")

    @property
    def set(self):
        return Interval(-oo, oo)

    def pdf(self, x):
        """
        Probability density function.

        Only available in closed form for special parameter values.
        For general parameters, use numerical methods.
        """
        alpha, beta, scale, location = self.alpha, self.beta, self.scale, self.location
        z = (x - location) / scale

        # Special case: alpha = 2 (Gaussian)
        # N(location, 2*scale^2)
        if alpha == 2:
            return exp(-z**2 / 2) / (scale * sqrt(2*pi))

        # Special case: alpha = 1, beta = 0 (Cauchy)
        if alpha == 1 and beta == 0:
            return 1 / (pi * scale * (1 + z**2))

        # Special case: alpha = 0.5, beta = 1 (Levy)
        if alpha == Rational(1, 2) and beta == 1:
            return Piecewise(
                (sqrt(scale / (2*pi)) * exp(-scale / (2*z)) / z**Rational(3, 2), z > 0),
                (0, True)
            )

        # No general closed form
        raise NotImplementedError(
            "No closed-form PDF for alpha=%s, beta=%s. "
            "Available special cases: alpha=2, (alpha=1, beta=0), (alpha=0.5, beta=1)." % (alpha, beta)
        )

    def _characteristic_function(self, t):
        r"""
        Characteristic function of the alpha-stable distribution.

        The characteristic function is defined as:

        For alpha != 1:
            phi(t) = exp(i*t*location - |scale*t|^alpha * (1 - i*beta*sign(t)*tan(pi*alpha/2)))

        For alpha = 1:
            phi(t) = exp(i*t*location - scale*|t| * (1 + i*beta*(2/pi)*sign(t)*log|t|))

        Parameters
        ==========

        t : Symbol
            The argument of the characteristic function.

        Returns
        =======

        Expr
            The characteristic function evaluated at t.
        """
        alpha, beta, scale, location = self.alpha, self.beta, self.scale, self.location

        # Sign function that handles t=0
        sign_t = Piecewise(
            (1, t > 0),
            (-1, t < 0),
            (0, True)
        )

        # Handle the alpha = 1 case separately (logarithmic term)
        if alpha == 1:
            return exp(
                I*t*location
                - scale * Abs(t) * (1 + I*beta * (2/pi) * sign_t * log(Abs(t)))
            )

        # For alpha != 1, use standard form
        return exp(
            I*t*location
            - (scale * Abs(t))**alpha * (1 - I*beta * sign_t * tan(pi*alpha/2))
        )


def AlphaStable(name, alpha, beta=0, scale=1, location=0):
    r"""
    Create an alpha-stable random variable.

    Parameters
    ==========

    name : Symbol or string
        Name of the random variable.
    alpha : Real number, 0 < alpha <= 2
        Stability parameter. Controls the tail heaviness of the distribution.
    beta : Real number, -1 <= beta <= 1
        Skewness parameter (default: 0 = symmetric).
    scale : Positive real number
        Scale parameter (default: 1).
    location : Real number
        Location parameter (default: 0).

    Returns
    =======

    RandomSymbol
        A continuous random variable with alpha-stable distribution.

    Examples
    ========

    >>> from sympy.stats import AlphaStable, density
    >>> from sympy import Symbol
    >>> x = Symbol('x', real=True)

    Create a Cauchy distribution (alpha=1, beta=0):

    >>> X = AlphaStable('X', 1, 0, 1, 0)
    >>> density(X)(x)
    1/(pi*(x**2 + 1))

    Create a standard normal distribution (alpha=2):

    >>> Y = AlphaStable('Y', 2, 0, 1, 0)
    >>> density(Y)(x)
    sqrt(2)*exp(-x**2/2)/(2*sqrt(pi))

    Create a general stable distribution with heavy tails:

    >>> Z = AlphaStable('Z', alpha=1.5, beta=0.5, scale=2, location=1)

    See Also
    ========

    sympy.stats.Normal
    sympy.stats.Cauchy

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Stable_distribution
    .. [2] Nolan, J. P. (2020). Univariate stable distributions.

    """
    return rv(name, AlphaStableDistribution, (alpha, beta, scale, location))
