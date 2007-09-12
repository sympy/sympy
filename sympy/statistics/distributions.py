from sympy.core import *
from sympy.functions import sqrt, exp, erf
import random


class Sample(tuple):
    """
    Sample([x1, x2, x3, ...]) represents a collection of samples.
    Sample parameters like mean, variance and stddev can be accessed as
    properties.
    """
    def __new__(cls, sample):
        s = tuple.__new__(cls, sample)
        s.mean = mean = sum(s) / len(s)
        s.variance = sum([(x-mean)**2 for x in s]) / Integer(len(s))
        s.stddev = sqrt(s.variance)
        return s

    @property
    def median(self):
        raise NotImplementedError

    def __repr__(self):
        return "Sample([" + ", ".join([str(x) for x in self]) + "])"

    __str__ = __repr__



class ContinuousProbability:
    """Base class for continuous probability distributions"""

    def probability(s, a, b):
        """Calculate the probability that a random number x generated
        from the distribution satisfies a <= x <= b """
        return s.cdf(b) - s.cdf(a)

    def random(s, n=None):
        """
        random() -- generate a random number from the distribution.
        random(n) -- generate a Sample of n random numbers.
        """
        if n is None:
            return s._random()
        else:
            return Sample([s._random() for i in xrange(n)])


class Normal(ContinuousProbability):
    """
    Normal(mu, sigma) represents the normal or Gaussian distribution
    with mean value mu and standard deviation sigma.

    Example usage:

        >>> N = Normal(1, 2)
        >>> N.mean
        1
        >>> N.variance
        4
        >>> N.probability(-oo, 1)   # probability on an interval
        1/2
        >>> N.probability(1, oo)
        1/2
        >>> N.probability(-oo, oo)
        1
        >>> N.probability(-2, 2)
        erf((1/2)**(1/2))
        >>> _.evalf()
        0.682689492137086

    """
    def __init__(self, mu, sigma):
        self.mu = Basic.sympify(mu)
        self.sigma = Basic.sympify(sigma)

    def __repr__(self):
        return "Normal(%s, %s)" % (self.mu, self.sigma)

    __str__ = __repr__

    mean = property(lambda s: s.mu)
    median = property(lambda s: s.mu)
    mode = property(lambda s: s.mu)
    stddev = property(lambda s: s.sigma)
    variance = property(lambda s: s.sigma**2)

    def pdf(s, x):
        """Return the probability density function as an expression in x"""
        x = Basic.sympify(x)
        return 1/(s.sigma*sqrt(2*pi)) * exp(-(x-s.mu)**2 / (2*s.sigma**2))

    def cdf(s, x):
        """Return the cumulative density function as an expression in x"""
        x = Basic.sympify(x)
        return (1+erf((x-s.mu)/(s.sigma*sqrt(2))))/2

    def _random(s):
        return random.gauss(float(s.mu), float(s.sigma))

    def confidence(s, p):
        """Return a symmetric (p*100)% confidence interval. For example,
        p=0.95 gives a 95% confidence interval. Currently this function
        only handles numerical values except in the trivial case p=1.

        Examples usage:
            # One standard deviation
            >>> N = Normal(0, 1)
            >>> N.confidence(0.68)
            (-1.95996398454005, 1.95996398454005)
            >>> N.probability(*_).evalf()
            0.68

            # Two standard deviations
            >>> N = Normal(0, 1)
            >>> N.confidence(0.95)
            (-1.95996398454005, 1.95996398454005)
            >>> N.probability(*_).evalf()
            0.95
        """

        if p == 1:
            return (-oo, oo)

        assert p <= 1

        # In terms of n*sigma, we have n = sqrt(2)*ierf(p). The inverse
        # error function is not yet implemented in SymPy but can easily be
        # computed numerically

        from sympy.numerics import Float, secant, evalf
        from sympy.numerics.functions2 import erf
        p = evalf(p)
        # calculate y = ierf(p) by solving erf(y) - p = 0
        y = secant(lambda y: erf(y) - p, 0)
        t = Real(str(evalf(s.sigma) * Float(2)**0.5 * y))
        mu = s.mu.evalf()
        return (mu-t, mu+t)

    @staticmethod
    def fit(sample):
        """Create a normal distribution fit to the mean and standard
        deviation of the given distribution or sample."""
        if not hasattr(sample, "stddev"):
            sample = Sample(sample)
        return Normal(sample.mean, sample.stddev)


class Uniform(ContinuousProbability):
    """
    Uniform(a, b) represents a probability distribution with uniform
    probability density on the interval [a, b] and zero density
    everywhere else.
    """
    def __init__(self, a, b):
        self.a = Basic.sympify(a)
        self.b = Basic.sympify(b)

    def __repr__(self):
        return "Uniform(%s, %s)" % (self.a, self.b)

    __str__ = __repr__

    mean = property(lambda s: (s.a+s.b)/2)
    median = property(lambda s: (s.a+s.b)/2)
    mode = property(lambda s: (s.a+s.b)/2)  # arbitrary
    variance = property(lambda s: (s.b-s.a)**2 / 12)
    stddev = property(lambda s: sqrt(s.variance))

    def pdf(s, x):
        """Return the probability density function as an expression in x"""
        x = Basic.sympify(x)
        if not isinstance(x, Number):
            raise NotImplementedError("SymPy does not yet support"
                "piecewise functions")
        if x < s.a or x > s.b:
            return Rational(0)
        return 1/(s.b-s.a)

    def cdf(s, x):
        """Return the cumulative density function as an expression in x"""
        x = Basic.sympify(x)
        if not isinstance(x, Number):
            raise NotImplementedError("SymPy does not yet support"
                "piecewise functions")
        if x <= s.a:
            return Rational(0)
        if x >= s.b:
            return Rational(1)
        return (x-s.a)/(s.b-s.a)

    def _random(s):
        return Real(random.uniform(float(s.a), float(s.b)))

    def confidence(s, p):
        """Generate a symmetric (p*100)% confidence interval.

        >>> U = Uniform(1, 2)
        >>> U.confidence(1)
        (1, 2)
        >>> U.confidence(Rational(1,2))
        (5/4, 7/4)
        """
        p = Basic.sympify(p)
        assert p <= 1

        d = (s.b-s.a)*p / 2
        return (s.mean - d, s.mean + d)

    @staticmethod
    def fit(sample):
        """Create a uniform distribution fit to the mean and standard
        deviation of the given distribution or sample."""
        if not hasattr(sample, "stddev"):
            sample = Sample(sample)
        m = sample.mean
        d = sqrt(12*sample.variance)/2
        return Uniform(m-d, m+d)
