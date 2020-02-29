"""
Finite Discrete Random Variables - Prebuilt variable types

Contains
========
FiniteRV
DiscreteUniform
Die
Bernoulli
Coin
Binomial
BetaBinomial
Hypergeometric
Rademacher
"""

from __future__ import print_function, division

import random

from sympy import (S, sympify, Rational, binomial, cacheit, Integer,
                   Dummy, Eq, Intersection, Interval,
                   Symbol, Lambda, Piecewise, Or, Gt, Lt, Ge, Le, Contains)
from sympy import beta as beta_fn
from sympy.external import import_module
from sympy.tensor.array import ArrayComprehensionMap
from sympy.stats.frv import (SingleFiniteDistribution,
                             SingleFinitePSpace)
from sympy.stats.rv import _value_check, Density, RandomSymbol

numpy = import_module('numpy')
scipy = import_module('scipy')
pymc3 = import_module('pymc3')

__all__ = ['FiniteRV',
'DiscreteUniform',
'Die',
'Bernoulli',
'Coin',
'Binomial',
'BetaBinomial',
'Hypergeometric',
'Rademacher'
]

def rv(name, cls, *args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return SingleFinitePSpace(name, dist).value

class FiniteDistributionHandmade(SingleFiniteDistribution):

    @property
    def dict(self):
        return self.args[0]

    def pmf(self, x):
        x = Symbol('x')
        return Lambda(x, Piecewise(*(
            [(v, Eq(k, x)) for k, v in self.dict.items()] + [(S.Zero, True)])))

    @property
    def set(self):
        return set(self.dict.keys())

    @staticmethod
    def check(density):
        for p in density.values():
            _value_check((p >= 0, p <= 1),
                        "Probability at a point must be between 0 and 1.")
        _value_check(Eq(sum(density.values()), 1), "Total Probability must be 1.")

def FiniteRV(name, density):
    r"""
    Create a Finite Random Variable given a dict representing the density.

    Parameters
    ==========

    density: A dict
        Dictionary conatining the pdf of finite distribution

    Examples
    ========

    >>> from sympy.stats import FiniteRV, P, E

    >>> density = {0: .1, 1: .2, 2: .3, 3: .4}
    >>> X = FiniteRV('X', density)

    >>> E(X)
    2.00000000000000
    >>> P(X >= 2)
    0.700000000000000

    Returns
    =======

    RandomSymbol

    """
    return rv(name, FiniteDistributionHandmade, density)

class DiscreteUniformDistribution(SingleFiniteDistribution):

    @staticmethod
    def check(*args):
        # not using _value_check since there is a
        # suggestion for the user
        if len(set(args)) != len(args):
            from sympy.utilities.iterables import multiset
            from sympy.utilities.misc import filldedent
            weights = multiset(args)
            n = Integer(len(args))
            for k in weights:
                weights[k] /= n
            raise ValueError(filldedent("""
                Repeated args detected but set expected. For a
                distribution having different weights for each
                item use the following:""") + (
                '\nS("FiniteRV(%s, %s)")' % ("'X'", weights)))

    @property
    def p(self):
        return Rational(1, len(self.args))

    @property  # type: ignore
    @cacheit
    def dict(self):
        return dict((k, self.p) for k in self.set)

    @property
    def set(self):
        return set(self.args)

    def pmf(self, x):
        if x in self.args:
            return self.p
        else:
            return S.Zero

    def _sample_random(self, size):
        x = Symbol('x')
        return ArrayComprehensionMap(lambda: self.args[random.randint(0, len(self.args)-1)], (x, 0, size)).doit()



def DiscreteUniform(name, items):
    r"""
    Create a Finite Random Variable representing a uniform distribution over
    the input set.

    Parameters
    ==========

    items: list/tuple
        Items over which Uniform distribution is to be made

    Examples
    ========

    >>> from sympy.stats import DiscreteUniform, density
    >>> from sympy import symbols

    >>> X = DiscreteUniform('X', symbols('a b c')) # equally likely over a, b, c
    >>> density(X).dict
    {a: 1/3, b: 1/3, c: 1/3}

    >>> Y = DiscreteUniform('Y', list(range(5))) # distribution over a range
    >>> density(Y).dict
    {0: 1/5, 1: 1/5, 2: 1/5, 3: 1/5, 4: 1/5}

    Returns
    =======

    RandomSymbol

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Discrete_uniform_distribution
    .. [2] http://mathworld.wolfram.com/DiscreteUniformDistribution.html

    """
    return rv(name, DiscreteUniformDistribution, *items)


class DieDistribution(SingleFiniteDistribution):
    _argnames = ('sides',)

    @staticmethod
    def check(sides):
        _value_check((sides.is_positive, sides.is_integer),
                    "number of sides must be a positive integer.")

    @property
    def is_symbolic(self):
        return not self.sides.is_number

    @property
    def high(self):
        return self.sides

    @property
    def low(self):
        return S.One

    @property
    def set(self):
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(0, self.sides))
        return set(map(Integer, list(range(1, self.sides + 1))))

    def pmf(self, x):
        x = sympify(x)
        if not (x.is_number or x.is_Symbol or isinstance(x, RandomSymbol)):
            raise ValueError("'x' expected as an argument of type 'number' or 'Symbol' or , "
                        "'RandomSymbol' not %s" % (type(x)))
        cond = Ge(x, 1) & Le(x, self.sides) & Contains(x, S.Integers)
        return Piecewise((S.One/self.sides, cond), (S.Zero, True))

def Die(name, sides=6):
    r"""
    Create a Finite Random Variable representing a fair die.

    Parameters
    ==========

    sides: Integer
        Represents the number of sides of the Die, by default is 6

    Examples
    ========

    >>> from sympy.stats import Die, density
    >>> from sympy import Symbol

    >>> D6 = Die('D6', 6) # Six sided Die
    >>> density(D6).dict
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}

    >>> D4 = Die('D4', 4) # Four sided Die
    >>> density(D4).dict
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}

    >>> n = Symbol('n', positive=True, integer=True)
    >>> Dn = Die('Dn', n) # n sided Die
    >>> density(Dn).dict
    Density(DieDistribution(n))
    >>> density(Dn).dict.subs(n, 4).doit()
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}

    Returns
    =======

    RandomSymbol
    """

    return rv(name, DieDistribution, sides)


class BernoulliDistribution(SingleFiniteDistribution):
    _argnames = ('p', 'succ', 'fail')

    @staticmethod
    def check(p, succ, fail):
        _value_check((p >= 0, p <= 1),
                    "p should be in range [0, 1].")

    @property
    def set(self):
        return set([self.succ, self.fail])

    def pmf(self, x):
        return Piecewise((self.p, x == self.succ),
                         (1 - self.p, x == self.fail),
                         (S.Zero, True))


def Bernoulli(name, p, succ=1, fail=0):
    r"""
    Create a Finite Random Variable representing a Bernoulli process.

    Parameters
    ==========

    p : Rational number between 0 and 1
       Represents probability of success
    succ : Integer/symbol/string
       Represents event of success
    fail : Integer/symbol/string
       Represents event of failure

    Examples
    ========

    >>> from sympy.stats import Bernoulli, density
    >>> from sympy import S

    >>> X = Bernoulli('X', S(3)/4) # 1-0 Bernoulli variable, probability = 3/4
    >>> density(X).dict
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli('X', S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> density(X).dict
    {Heads: 1/2, Tails: 1/2}

    Returns
    =======

    RandomSymbol

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Bernoulli_distribution
    .. [2] http://mathworld.wolfram.com/BernoulliDistribution.html

    """

    return rv(name, BernoulliDistribution, p, succ, fail)


def Coin(name, p=S.Half):
    r"""
    Create a Finite Random Variable representing a Coin toss.

    Parameters
    ==========

    p : Rational Numeber between 0 and 1
      Represents probability of getting "Heads", by default is Half

    Examples
    ========

    >>> from sympy.stats import Coin, density
    >>> from sympy import Rational

    >>> C = Coin('C') # A fair coin toss
    >>> density(C).dict
    {H: 1/2, T: 1/2}

    >>> C2 = Coin('C2', Rational(3, 5)) # An unfair coin
    >>> density(C2).dict
    {H: 3/5, T: 2/5}

    Returns
    =======

    RandomSymbol

    See Also
    ========

    sympy.stats.Binomial

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Coin_flipping

    """
    return rv(name, BernoulliDistribution, p, 'H', 'T')


class BinomialDistribution(SingleFiniteDistribution):
    _argnames = ('n', 'p', 'succ', 'fail')

    @staticmethod
    def check(n, p, succ, fail):
        _value_check((n.is_integer, n.is_nonnegative),
                    "'n' must be nonnegative integer.")
        _value_check((p <= 1, p >= 0),
                    "p should be in range [0, 1].")

    @property
    def high(self):
        return self.n

    @property
    def low(self):
        return S.Zero

    @property
    def is_symbolic(self):
        return not self.n.is_number

    @property
    def set(self):
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(0, self.n))
        return set(self.dict.keys())

    def pmf(self, x):
        n, p = self.n, self.p
        x = sympify(x)
        if not (x.is_number or x.is_Symbol or isinstance(x, RandomSymbol)):
            raise ValueError("'x' expected as an argument of type 'number' or 'Symbol' or , "
                        "'RandomSymbol' not %s" % (type(x)))
        cond = Ge(x, 0) & Le(x, n) & Contains(x, S.Integers)
        return Piecewise((binomial(n, x) * p**x * (1 - p)**(n - x), cond), (S.Zero, True))

    @property  # type: ignore
    @cacheit
    def dict(self):
        if self.is_symbolic:
            return Density(self)
        return dict((k*self.succ + (self.n-k)*self.fail, self.pmf(k))
                    for k in range(0, self.n + 1))

def Binomial(name, n, p, succ=1, fail=0):
    r"""
    Create a Finite Random Variable representing a binomial distribution.

    Parameters
    ==========

    n : Positive Integer
      Represents number of trials
    p : Rational Number between 0 and 1
      Represents probability of success
    succ : Integer/symbol/string
      Represents event of success, by default is 1
    fail : Integer/symbol/string
      Represents event of failure, by default is 0

    Examples
    ========

    >>> from sympy.stats import Binomial, density
    >>> from sympy import S, Symbol

    >>> X = Binomial('X', 4, S.Half) # Four "coin flips"
    >>> density(X).dict
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}

    >>> n = Symbol('n', positive=True, integer=True)
    >>> p = Symbol('p', positive=True)
    >>> X = Binomial('X', n, S.Half) # n "coin flips"
    >>> density(X).dict
    Density(BinomialDistribution(n, 1/2, 1, 0))
    >>> density(X).dict.subs(n, 4).doit()
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}

    Returns
    =======

    RandomSymbol

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Binomial_distribution
    .. [2] http://mathworld.wolfram.com/BinomialDistribution.html

    """

    return rv(name, BinomialDistribution, n, p, succ, fail)

#-------------------------------------------------------------------------------
# Beta-binomial distribution ----------------------------------------------------------

class BetaBinomialDistribution(SingleFiniteDistribution):
    _argnames = ('n', 'alpha', 'beta')

    @staticmethod
    def check(n, alpha, beta):
        _value_check((n.is_integer, n.is_nonnegative),
        "'n' must be nonnegative integer. n = %s." % str(n))
        _value_check((alpha > 0),
        "'alpha' must be: alpha > 0 . alpha = %s" % str(alpha))
        _value_check((beta > 0),
        "'beta' must be: beta > 0 . beta = %s" % str(beta))

    @property
    def high(self):
        return self.n

    @property
    def low(self):
        return S.Zero

    @property
    def is_symbolic(self):
        return not self.n.is_number

    @property
    def set(self):
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(0, self.n))
        return set(map(Integer, list(range(0, self.n + 1))))

    def pmf(self, k):
        n, a, b = self.n, self.alpha, self.beta
        return binomial(n, k) * beta_fn(k + a, n - k + b) / beta_fn(a, b)

    def _sample_pymc3(self, size):
        n, a, b = int(self.n), float(self.alpha), float(self.beta)
        with pymc3.Model():
            pymc3.BetaBinomial('X', alpha=a, beta=b, n=n)
            return pymc3.sample(size, chains=1, progressbar=False)[:]['X']

def BetaBinomial(name, n, alpha, beta):
    r"""
    Create a Finite Random Variable representing a Beta-binomial distribution.

    Parameters
    ==========

    n : Positive Integer
      Represents number of trials
    alpha : Real positive number
    beta : Real positive number

    Examples
    ========

    >>> from sympy.stats import BetaBinomial, density

    >>> X = BetaBinomial('X', 2, 1, 1)
    >>> density(X).dict
    {0: 1/3, 1: 2*beta(2, 2), 2: 1/3}

    Returns
    =======

    RandomSymbol

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Beta-binomial_distribution
    .. [2] http://mathworld.wolfram.com/BetaBinomialDistribution.html

    """

    return rv(name, BetaBinomialDistribution, n, alpha, beta)


class HypergeometricDistribution(SingleFiniteDistribution):
    _argnames = ('N', 'm', 'n')

    @staticmethod
    def check(n, N, m):
        _value_check((N.is_integer, N.is_nonnegative),
                     "'N' must be nonnegative integer. N = %s." % str(n))
        _value_check((n.is_integer, n.is_nonnegative),
                     "'n' must be nonnegative integer. n = %s." % str(n))
        _value_check((m.is_integer, m.is_nonnegative),
                     "'m' must be nonnegative integer. m = %s." % str(n))

    @property
    def is_symbolic(self):
        return any(not x.is_number for x in (self.N, self.m, self.n))

    @property
    def high(self):
        return Piecewise((self.n, Lt(self.n, self.m) != False), (self.m, True))

    @property
    def low(self):
        return Piecewise((0, Gt(0, self.n + self.m - self.N) != False), (self.n + self.m - self.N, True))

    @property
    def set(self):
        N, m, n = self.N, self.m, self.n
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(self.low, self.high))
        return set([i for i in range(max(0, n + m - N), min(n, m) + 1)])

    def pmf(self, k):
        N, m, n = self.N, self.m, self.n
        return S(binomial(m, k) * binomial(N - m, n - k))/binomial(N, n)

    def _sample_scipy(self, size):
        import scipy.stats # Make sure that stats is imported
        N, m, n = int(self.N), int(self.m), int(self.n)
        return scipy.stats.hypergeom.rvs(M=m, n=n, N=N, size=size)

def Hypergeometric(name, N, m, n):
    r"""
    Create a Finite Random Variable representing a hypergeometric distribution.

    Parameters
    ==========

    N : Positive Integer
      Represents finite population of size N.
    m : Positive Integer
      Represents number of trials with required feature.
    n : Positive Integer
      Represents numbers of draws.


    Examples
    ========

    >>> from sympy.stats import Hypergeometric, density

    >>> X = Hypergeometric('X', 10, 5, 3) # 10 marbles, 5 white (success), 3 draws
    >>> density(X).dict
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}

    Returns
    =======

    RandomSymbol

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hypergeometric_distribution
    .. [2] http://mathworld.wolfram.com/HypergeometricDistribution.html

    """
    return rv(name, HypergeometricDistribution, N, m, n)


class RademacherDistribution(SingleFiniteDistribution):

    @property
    def set(self):
        return set([-1, 1])

    @property
    def pmf(self):
        k = Dummy('k')
        return Lambda(k, Piecewise((S.Half, Or(Eq(k, -1), Eq(k, 1))), (S.Zero, True)))

def Rademacher(name):
    r"""
    Create a Finite Random Variable representing a Rademacher distribution.

    Examples
    ========

    >>> from sympy.stats import Rademacher, density

    >>> X = Rademacher('X')
    >>> density(X).dict
    {-1: 1/2, 1: 1/2}

    Returns
    =======

    RandomSymbol

    See Also
    ========

    sympy.stats.Bernoulli

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Rademacher_distribution

    """
    return rv(name, RademacherDistribution)
