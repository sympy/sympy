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

from sympy import (S, sympify, Rational, binomial, cacheit, Integer,
        Dict, Basic, KroneckerDelta, Dummy, Eq, Intersection, Interval,
        Symbol, Lambda, Piecewise, Or)
from sympy import beta as beta_fn
from sympy.concrete.summations import Sum
from sympy.core.compatibility import as_int, range
from sympy.stats.rv import _value_check, Density
from sympy.stats.frv import (SingleFinitePSpace, SingleFiniteDistribution)

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

def _is_sym_dim(n):
    n = sympify(n)
    if n.atoms(Symbol, Integer) == n.atoms() and\
        n.has(Symbol):
        return True
    return False

def rv(name, cls, *args):
    args = list(map(sympify, args))
    i = 0
    while i < len(args): # Converting to Dict since dict is not hashable
        if isinstance(args[i], dict):
            args[i] = Dict(args[i])
        i += 1
    dist = cls(*args)
    dist.check(*args)
    return SingleFinitePSpace(name, dist).value

class FiniteDistributionHandmade(SingleFiniteDistribution):

    @property
    def dict(self):
        return self.args[0]

    def pdf(self, x):
        x = Symbol('x')
        return Lambda(x, Piecewise(*(
            [(v, Eq(k, x)) for k, v in self.dict.items()] + [(0, True)])))

    @property
    def set(self):
        return list(self.dict.keys())

    @staticmethod
    def check(density):
        for p in density.values():
            _value_check((p >= 0, p <= 1),
                        "Probability at a point must be between 0 and 1.")
        _value_check(Eq(sum(density.values()), 1), "Total Probability must be 1.")

def FiniteRV(name, density):
    """
    Create a Finite Random Variable given a dict representing the density.

    Returns a RandomSymbol.

    >>> from sympy.stats import FiniteRV, P, E

    >>> density = {0: .1, 1: .2, 2: .3, 3: .4}
    >>> X = FiniteRV('X', density)

    >>> E(X)
    2.00000000000000
    >>> P(X >= 2)
    0.700000000000000
    """
    return rv(name, FiniteDistributionHandmade, density)

class DiscreteUniformDistribution(SingleFiniteDistribution):
    @property
    def p(self):
        return Rational(1, len(self.args))

    @property
    @cacheit
    def dict(self):
        return dict((k, self.p) for k in self.set)

    @property
    def set(self):
        return list(self.args)

    def pdf(self, x):
        if x in self.args:
            return self.p
        else:
            return S.Zero


def DiscreteUniform(name, items):
    """
    Create a Finite Random Variable representing a uniform distribution over
    the input set.

    Returns a RandomSymbol.

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
        return _is_sym_dim(self.sides)

    @property
    def high(self):
        return self.sides

    @property
    def low(self):
        return S(0)

    @property
    def set(self):
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(0, self.sides))
        return list(map(Integer, list(range(1, self.sides + 1))))

    def pdf(self, x):
        x = sympify(x)
        if x.is_number:
            if x.is_Integer and x >= 1 and x <= self.sides:
                return Rational(1, self.sides)
            return S.Zero
        if x.is_Symbol:
            i = Dummy('i', integer=True, positive=True)
            t = Dummy('t')
            return Sum(KroneckerDelta(x, i)/self.sides, (i, 1, self.sides))
        raise ValueError("'x' expected as an argument of type 'number' or 'symbol', "
                        "not %s" % (type(x)))


def Die(name, sides=6):
    """
    Create a Finite Random Variable representing a fair die.

    Returns a RandomSymbol.

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
        return [self.succ, self.fail]

    def pdf(self, x):
        return Piecewise((self.p, x == self.succ), (1 - self.p, x == self.fail), (0, True))


def Bernoulli(name, p, succ=1, fail=0):
    """
    Create a Finite Random Variable representing a Bernoulli process.

    Returns a RandomSymbol

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

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Bernoulli_distribution
    .. [2] http://mathworld.wolfram.com/BernoulliDistribution.html

    """

    return rv(name, BernoulliDistribution, p, succ, fail)


def Coin(name, p=S.Half):
    """
    Create a Finite Random Variable representing a Coin toss.

    Probability p is the chance of gettings "Heads." Half by default

    Returns a RandomSymbol.

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
        return S(0)

    @property
    def is_symbolic(self):
        return _is_sym_dim(self.n)

    @property
    def set(self):
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(0, self.n))
        return list(self.dict.keys())


    def pdf(self, k):
        n, p = self.n, self.p
        return binomial(n, k) * p**k * (1 - p)**(n - k)

    @property
    @cacheit
    def dict(self):
        if self.is_symbolic:
            return Density(self)
        return dict((k*self.succ + (self.n-k)*self.fail, self.pdf(k))
                    for k in range(0, self.n + 1))

def Binomial(name, n, p, succ=1, fail=0):
    """
    Create a Finite Random Variable representing a binomial distribution.

    Returns a RandomSymbol.

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
        return S(0)

    @property
    def is_symbolic(self):
        return _is_sym_dim(self.n)

    @property
    def set(self):
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(0, self.n))
        return list(map(Integer, list(range(0, self.n + 1))))

    def pdf(self, k):
        n, a, b = self.n, self.alpha, self.beta
        return binomial(n, k) * beta_fn(k + a, n - k + b) / beta_fn(a, b)


def BetaBinomial(name, n, alpha, beta):
    """
    Create a Finite Random Variable representing a Beta-binomial distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import BetaBinomial, density
    >>> from sympy import S

    >>> X = BetaBinomial('X', 2, 1, 1)
    >>> density(X).dict
    {0: beta(1, 3)/beta(1, 1), 1: 2*beta(2, 2)/beta(1, 1), 2: beta(3, 1)/beta(1, 1)}

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Beta-binomial_distribution
    .. [2] http://mathworld.wolfram.com/BetaBinomialDistribution.html

    """

    return rv(name, BetaBinomialDistribution, n, alpha, beta)


class HypergeometricDistribution(SingleFiniteDistribution):
    _argnames = ('N', 'm', 'n')

    @property
    def is_symbolic(self):
        return any(_is_sym_dim(x) for x in (self.N, self.m, self.n))

    @property
    def high(self):
        return max(self.n, self.K)

    @property
    def low(self):
        return min(0, self.n + self.m - self.N)

    @property
    def set(self):
        N, m, n = self.N, self.m, self.n
        if self.is_symbolic:
            return Intersection(S.Naturals0, Interval(max(0, n + m - N), min(n, m)))
        return [i for i in range(max(0, n + m - N), min(n, m) + 1)]

    def pdf(self, k):
        N, m, n = self.N, self.m, self.n
        return Rational(binomial(m, k) * binomial(N - m, n - k), binomial(N, n))

def Hypergeometric(name, N, m, n):
    """
    Create a Finite Random Variable representing a hypergeometric distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Hypergeometric, density
    >>> from sympy import S

    >>> X = Hypergeometric('X', 10, 5, 3) # 10 marbles, 5 white (success), 3 draws
    >>> density(X).dict
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hypergeometric_distribution
    .. [2] http://mathworld.wolfram.com/HypergeometricDistribution.html

    """
    return rv(name, HypergeometricDistribution, N, m, n)


class RademacherDistribution(SingleFiniteDistribution):

    @property
    def set(self):
        return [-1, 1]

    @property
    def pdf(self):
        k = Dummy('k')
        return Lambda(k, Piecewise((S.Half, Or(Eq(k, -1), Eq(k, 1))), (0, True)))

def Rademacher(name):
    """
    Create a Finite Random Variable representing a Rademacher distribution.

    Return a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Rademacher, density

    >>> X = Rademacher('X')
    >>> density(X).dict
    {-1: 1/2, 1: 1/2}

    See Also
    ========

    sympy.stats.Bernoulli

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Rademacher_distribution

    """
    return rv(name, RademacherDistribution)
