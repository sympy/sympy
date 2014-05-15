from __future__ import print_function, division

from sympy.core import sympify, Lambda, Dummy, Integer, Rational, oo, Float, pi
from sympy.core.compatibility import xrange
from sympy.functions import sqrt, exp, erf
from sympy.printing import sstr
from sympy.utilities import default_sort_key

import random


class Sample(tuple):
    """
    Sample([x1, x2, x3, ...]) represents a collection of samples.
    Sample parameters like mean, variance and stddev can be accessed as
    properties.
    The sample will be sorted.

    Examples
    ========

        >>> from sympy.statistics.distributions import Sample
        >>> Sample([0, 1, 2, 3])
        Sample([0, 1, 2, 3])
        >>> Sample([8, 3, 2, 4, 1, 6, 9, 2])
        Sample([1, 2, 2, 3, 4, 6, 8, 9])
        >>> s = Sample([1, 2, 3, 4, 5])
        >>> s.mean
        3
        >>> s.stddev
        sqrt(2)
        >>> s.median
        3
        >>> s.variance
        2

    """
    def __new__(cls, sample):
        s = tuple.__new__(cls, sorted(sample, key=default_sort_key))
        s.mean = mean = sum(s) / Integer(len(s))
        s.variance = sum([(x - mean)**2 for x in s]) / Integer(len(s))
        s.stddev = sqrt(s.variance)
        if len(s) % 2:
            s.median = s[len(s)//2]
        else:
            s.median = sum(s[len(s)//2 - 1:len(s)//2 + 1]) / Integer(2)
        return s

    def __repr__(self):
        return sstr(self)

    def __str__(self):
        return sstr(self)


class ContinuousProbability(object):
    """Base class for continuous probability distributions"""

    def probability(s, a, b):
        """
        Calculate the probability that a random number x generated
        from the distribution satisfies a <= x <= b

        Examples
        ========

            >>> from sympy.statistics import Normal
            >>> from sympy.core import oo
            >>> Normal(0, 1).probability(-1, 1)
            erf(sqrt(2)/2)
            >>> Normal(0, 1).probability(1, oo)
            -erf(sqrt(2)/2)/2 + 1/2

        """
        return s.cdf(b) - s.cdf(a)

    def random(s, n=None):
        """
        random() -- generate a random number from the distribution.
        random(n) -- generate a Sample of n random numbers.

        Examples
        ========

            >>> from sympy.statistics import Uniform
            >>> x = Uniform(1, 5).random()
            >>> x < 5 and x > 1
            True
            >>> x = Uniform(-4, 2).random()
            >>> x < 2 and x > -4
            True

        """
        if n is None:
            return s._random()
        else:
            return Sample([s._random() for i in xrange(n)])

    def __repr__(self):
        return sstr(self)

    def __str__(self):
        return sstr(self)


class Normal(ContinuousProbability):
    """
    Normal(mu, sigma) represents the normal or Gaussian distribution
    with mean value mu and standard deviation sigma.

    Examples
    ========

        >>> from sympy.statistics import Normal
        >>> from sympy import oo
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
        >>> N.probability(-1, 3)
        erf(sqrt(2)/2)
        >>> _.evalf()
        0.682689492137086

    """
    def __init__(self, mu, sigma):
        self.mu = sympify(mu)
        self.sigma = sympify(sigma)

    mean = property(lambda s: s.mu)
    median = property(lambda s: s.mu)
    mode = property(lambda s: s.mu)
    stddev = property(lambda s: s.sigma)
    variance = property(lambda s: s.sigma**2)

    def pdf(s, x):
        """
        Return the probability density function as an expression in x

        Examples
        ========

            >>> from sympy.statistics import Normal
            >>> Normal(1, 2).pdf(0)
            sqrt(2)*exp(-1/8)/(4*sqrt(pi))
            >>> from sympy.abc import x
            >>> Normal(1, 2).pdf(x)
            sqrt(2)*exp(-(x - 1)**2/8)/(4*sqrt(pi))

        """
        x = sympify(x)
        return 1/(s.sigma*sqrt(2*pi)) * exp(-(x - s.mu)**2 / (2*s.sigma**2))

    def cdf(s, x):
        """
        Return the cumulative density function as an expression in x

        Examples
        ========

            >>> from sympy.statistics import Normal
            >>> Normal(1, 2).cdf(0)
            -erf(sqrt(2)/4)/2 + 1/2
            >>> from sympy.abc import x
            >>> Normal(1, 2).cdf(x)
            erf(sqrt(2)*(x - 1)/4)/2 + 1/2

        """
        x = sympify(x)
        return (1 + erf((x - s.mu)/(s.sigma*sqrt(2))))/2

    def _random(s):
        return random.gauss(float(s.mu), float(s.sigma))

    def confidence(s, p):
        """Return a symmetric (p*100)% confidence interval. For example,
        p=0.95 gives a 95% confidence interval. Currently this function
        only handles numerical values except in the trivial case p=1.

        For example, one standard deviation:

            >>> from sympy.statistics import Normal
            >>> N = Normal(0, 1)
            >>> N.confidence(0.68)
            (-0.994457883209753, 0.994457883209753)
            >>> N.probability(*_).evalf()
            0.680000000000000

        Two standard deviations:

            >>> N = Normal(0, 1)
            >>> N.confidence(0.95)
            (-1.95996398454005, 1.95996398454005)
            >>> N.probability(*_).evalf()
            0.950000000000000

        """

        if p == 1:
            return (-oo, oo)

        if p > 1:
            raise ValueError("p cannot be greater than 1")

        # In terms of n*sigma, we have n = sqrt(2)*ierf(p). The inverse
        # error function is not yet implemented in SymPy but can easily be
        # computed numerically

        from sympy.mpmath import mpf, erfinv

        # calculate y = ierf(p) by solving erf(y) - p = 0
        y = erfinv(mpf(p))
        t = Float(str(mpf(float(s.sigma)) * mpf(2)**0.5 * y))
        mu = s.mu.evalf()
        return (mu - t, mu + t)

    @staticmethod
    def fit(sample):
        """
        Create a normal distribution fit to the mean and standard
        deviation of the given distribution or sample.

        Examples
        ========

            >>> from sympy.statistics import Normal
            >>> Normal.fit([1,2,3,4,5])
            Normal(3, sqrt(2))
            >>> from sympy.abc import x, y
            >>> Normal.fit([x, y])
            Normal(x/2 + y/2, sqrt((-x/2 + y/2)**2/2 + (x/2 - y/2)**2/2))

        """
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
        self.a = sympify(a)
        self.b = sympify(b)

    mean = property(lambda s: (s.a + s.b)/2)
    median = property(lambda s: (s.a + s.b)/2)
    mode = property(lambda s: (s.a + s.b)/2)  # arbitrary
    variance = property(lambda s: (s.b - s.a)**2 / 12)
    stddev = property(lambda s: sqrt(s.variance))

    def pdf(s, x):
        """
        Return the probability density function as an expression in x

        Examples
        ========

            >>> from sympy.statistics import Uniform
            >>> Uniform(1, 5).pdf(1)
            1/4
            >>> Uniform(2, 4).pdf(2)
            1/2

        """
        x = sympify(x)
        if not x.is_Number:
            raise NotImplementedError("SymPy does not yet support "
                "piecewise functions")
        if x < s.a or x > s.b:
            return Rational(0)
        return 1/(s.b - s.a)

    def cdf(s, x):
        """
        Return the cumulative density function as an expression in x

        Examples
        ========

            >>> from sympy.statistics import Uniform
            >>> Uniform(1, 5).cdf(2)
            1/4
            >>> Uniform(1, 5).cdf(4)
            3/4

        """
        x = sympify(x)
        if not x.is_Number:
            raise NotImplementedError("SymPy does not yet support "
                "piecewise functions")
        if x <= s.a:
            return Rational(0)
        if x >= s.b:
            return Rational(1)
        return (x - s.a)/(s.b - s.a)

    def _random(s):
        return Float(random.uniform(float(s.a), float(s.b)))

    def confidence(s, p):
        """Generate a symmetric (p*100)% confidence interval.

        >>> from sympy import Rational
        >>> from sympy.statistics import Uniform
        >>> U = Uniform(1, 2)
        >>> U.confidence(1)
        (1, 2)
        >>> U.confidence(Rational(1,2))
        (5/4, 7/4)

        """
        p = sympify(p)
        if p > 1:
            raise ValueError("p cannot be greater than 1")

        d = (s.b - s.a)*p / 2
        return (s.mean - d, s.mean + d)

    @staticmethod
    def fit(sample):
        """
        Create a uniform distribution fit to the mean and standard
        deviation of the given distribution or sample.

        Examples
        ========

            >>> from sympy.statistics import Uniform
            >>> Uniform.fit([1, 2, 3, 4, 5])
            Uniform(-sqrt(6) + 3, sqrt(6) + 3)
            >>> Uniform.fit([1, 2])
            Uniform(-sqrt(3)/2 + 3/2, sqrt(3)/2 + 3/2)

        """
        if not hasattr(sample, "stddev"):
            sample = Sample(sample)
        m = sample.mean
        d = sqrt(12*sample.variance)/2
        return Uniform(m - d, m + d)


class PDF(ContinuousProbability):
    """
    PDF(func, (x, a, b)) represents continuous probability distribution
    with probability distribution function func(x) on interval (a, b)

    If func is not normalized so that integrate(func, (x, a, b)) == 1,
    it can be normalized using PDF.normalize() method

    Examples
    ========

        >>> from sympy import Symbol, exp, oo
        >>> from sympy.statistics.distributions import PDF
        >>> from sympy.abc import x
        >>> a = Symbol('a', positive=True)

        >>> exponential = PDF(exp(-x/a)/a, (x,0,oo))
        >>> exponential.pdf(x)
        exp(-x/a)/a
        >>> exponential.cdf(x)
        1 - exp(-x/a)
        >>> exponential.mean
        a
        >>> exponential.variance
        a**2

    """

    def __init__(self, pdf, var):
        #XXX maybe add some checking of parameters
        if isinstance(var, (tuple, list)):
            self.pdf = Lambda(var[0], pdf)
            self.domain = tuple(var[1:])
        else:
            self.pdf = Lambda(var, pdf)
            self.domain = (-oo, oo)
        self._cdf = None
        self._mean = None
        self._variance = None
        self._stddev = None

    def normalize(self):
        """
        Normalize the probability distribution function so that
        integrate(self.pdf(x), (x, a, b)) == 1

        Examples
        ========

            >>> from sympy import Symbol, exp, oo
            >>> from sympy.statistics.distributions import PDF
            >>> from sympy.abc import x
            >>> a = Symbol('a', positive=True)

            >>> exponential = PDF(exp(-x/a), (x,0,oo))
            >>> exponential.normalize().pdf(x)
            exp(-x/a)/a

        """

        norm = self.probability(*self.domain)
        if norm != 1:
            w = Dummy('w', real=True)
            return self.__class__(self.pdf(w)/norm, (w, self.domain[0], self.domain[1]))
            #self._cdf = Lambda(w, (self.cdf(w) - self.cdf(self.domain[0]))/norm)
            #if self._mean is not None:
            #    self._mean /= norm
            #if self._variance is not None:
            #    self._variance = (self._variance + (self._mean*norm)**2)/norm - self.mean**2
            #if self._stddev is not None:
            #    self._stddev = sqrt(self._variance)
        else:
            return self

    def cdf(self, x):
        """
        Return the cumulative density function as an expression in x

        Examples
        ========

            >>> from sympy.statistics.distributions import PDF
            >>> from sympy import exp, oo
            >>> from sympy.abc import x, y
            >>> PDF(exp(-x/y), (x,0,oo)).cdf(4)
            y - y*exp(-4/y)
            >>> PDF(2*x + y, (x, 10, oo)).cdf(0)
            -10*y - 100

        """
        x = sympify(x)
        if self._cdf is not None:
            return self._cdf(x)
        else:
            from sympy import integrate
            w = Dummy('w', real=True)
            self._cdf = integrate(self.pdf(w), w)
            self._cdf = Lambda(
                w, self._cdf - self._cdf.subs(w, self.domain[0]))
            return self._cdf(x)

    def _get_mean(self):
        if self._mean is not None:
            return self._mean
        else:
            from sympy import integrate
            w = Dummy('w', real=True)
            self._mean = integrate(
                self.pdf(w)*w, (w, self.domain[0], self.domain[1]))
            return self._mean

    def _get_variance(self):
        if self._variance is not None:
            return self._variance
        else:
            from sympy import integrate, simplify
            w = Dummy('w', real=True)
            self._variance = integrate(self.pdf(
                w)*w**2, (w, self.domain[0], self.domain[1])) - self.mean**2
            self._variance = simplify(self._variance)
            return self._variance

    def _get_stddev(self):
        if self._stddev is not None:
            return self._stddev
        else:
            self._stddev = sqrt(self.variance)
            return self._stddev

    mean = property(_get_mean)
    variance = property(_get_variance)
    stddev = property(_get_stddev)

    def _random(s):
        raise NotImplementedError

    def transform(self, func, var):
        """
        Return a probability distribution of random variable func(x)
        currently only some simple injective functions are supported

        Examples
        ========

            >>> from sympy.statistics.distributions import PDF
            >>> from sympy import oo
            >>> from sympy.abc import x, y
            >>> PDF(2*x + y, (x, 10, oo)).transform(x, y)
            PDF(0, ((_w,), x, x))

        """

        w = Dummy('w', real=True)

        from sympy import solve
        from sympy import S
        inverse = solve(func - w, var)
        newPdf = S.Zero
        funcdiff = func.diff(var)
        #TODO check if x is in domain
        for x in inverse:
            # this assignment holds only for x in domain
            # in general it would require implementing
            # piecewise defined functions in sympy
            newPdf += (self.pdf(var)/abs(funcdiff)).subs(var, x)

        return PDF(newPdf, (w, func.subs(var, self.domain[0]), func.subs(var, self.domain[1])))
