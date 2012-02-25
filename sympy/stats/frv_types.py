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
Hypergeometric
"""

from sympy.stats.frv import SingleFinitePSpace, create_SingleFinitePSpace
from sympy import S, sympify, Rational, binomial

__all__ = ['FiniteRV', 'DiscreteUniform', 'Die', 'Bernoulli', 'Coin',
        'Binomial', 'Hypergeometric']

def FiniteRV(density, symbol=None):
    """
    Create a Finite Random Variable given a dict representing the density.

    Returns a RandomSymbol.

    >>> from sympy.stats import FiniteRV, P, E

    >>> density = {0: .1, 1: .2, 2: .3, 3: .4}
    >>> X = FiniteRV(density)

    >>> E(X)
    2.00000000000000
    >>> P(X>=2)
    0.700000000000000
    """
    return create_SingleFinitePSpace(density, symbol).value

class DiscreteUniformPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a discrete uniform
    distribution.

    This class is for internal use.

    Create DiscreteUniform Random Symbols using DiscreteUniform function

    Examples
    ========

    >>> from sympy.stats import DiscreteUniform, Density
    >>> from sympy import symbols

    >>> X = DiscreteUniform(symbols('a b c')) # equally likely over a, b, c
    >>> Density(X)
    {a: 1/3, b: 1/3, c: 1/3}

    >>> Y = DiscreteUniform(range(5)) # distribution over a range
    >>> Density(Y)
    {0: 1/5, 1: 1/5, 2: 1/5, 3: 1/5, 4: 1/5}
    """
    def __new__(cls, items, symbol=None):
        density = dict((sympify(item), Rational(1, len(items)))
                       for item in items)
        return create_SingleFinitePSpace(density, symbol, cls)

def DiscreteUniform(items, symbol=None):
    """
    Create a Finite Random Variable representing a uniform distribution over
    the input set.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import DiscreteUniform, Density
    >>> from sympy import symbols

    >>> X = DiscreteUniform(symbols('a b c')) # equally likely over a, b, c
    >>> Density(X)
    {a: 1/3, b: 1/3, c: 1/3}

    >>> Y = DiscreteUniform(range(5)) # distribution over a range
    >>> Density(Y)
    {0: 1/5, 1: 1/5, 2: 1/5, 3: 1/5, 4: 1/5}

    """
    return DiscreteUniformPSpace(items, symbol).value

class DiePSpace(DiscreteUniformPSpace):
    """
    Create a Finite Random Variable representing a fair die.

    This class is for internal use.

    Create Dice Random Symbols using Die function

    >>> from sympy.stats import Die, Density

    >>> X = Die(6) # Six sided Die
    >>> Density(X)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}

    >>> X = Die(4) # Four sided Die
    >>> Density(X)
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}
    """
    _count = 0
    _name = 'die'
    def __new__(cls, sides=6, symbol=None):
        return DiscreteUniformPSpace.__new__(cls, range(1, sides+1), symbol)

def Die(sides=6, symbol=None):
    """
    Create a Finite Random Variable representing a fair die.

    Returns a RandomSymbol.

    >>> from sympy.stats import Die, Density

    >>> X = Die(6) # Six sided Die
    >>> Density(X)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}

    >>> X = Die(4) # Four sided Die
    >>> Density(X)
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}
    """

    return DiePSpace(sides, symbol).value

class BernoulliPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a Bernoulli process.

    Returns a RandomSymbol.

    This class is for internal use.

    Create Bernoulli Random Symbols using Bernoulli function.

    >>> from sympy.stats import Bernoulli, Density
    >>> from sympy import S

    >>> X = Bernoulli(S(3)/4) # 1-0 Bernoulli variable, probability = 3/4
    >>> Density(X)
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli(S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> Density(X)
    {Heads: 1/2, Tails: 1/2}
    """

    _name = 'bernoulli'
    def __new__(cls, p, succ=1, fail=0, symbol=None):
        succ, fail, p = map(sympify, (succ, fail, p))
        density = {succ:p, fail:(1-p)}
        return create_SingleFinitePSpace(density, symbol, cls)

def Bernoulli(p, succ=1, fail=0, symbol=None):
    """
    Create a Finite Random Variable representing a Bernoulli process.

    Returns a RandomSymbol

    >>> from sympy.stats import Bernoulli, Density
    >>> from sympy import S

    >>> X = Bernoulli(S(3)/4, 1, 0) # 1-0 Bernoulli variable, probability = 3/4
    >>> Density(X)
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli(S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> Density(X)
    {Heads: 1/2, Tails: 1/2}
    """

    return BernoulliPSpace(p, succ, fail, symbol).value

class CoinPSpace(BernoulliPSpace):
    """
    A probability space representing a coin toss.

    Probability p is the chance of gettings "Heads." Half by default

    This class is for internal use.

    Create Coin's using Coin function

    >>> from sympy.stats import Coin, Density
    >>> from sympy import Rational

    >>> X = Coin() # A fair coin toss
    >>> Density(X)
    {H: 1/2, T: 1/2}

    >>> X = Coin(Rational(3, 5)) # An unfair coin
    >>> Density(X)
    {H: 3/5, T: 2/5}
    """
    _count = 0
    _name = 'coin'
    def __new__(cls, p=S.Half, symbol=None):
        return BernoulliPSpace.__new__(cls, p, 'H', 'T', symbol)

def Coin(p=S.Half, symbol=None):
    """
    Create a Finite Random Variable representing a Coin toss.

    Probability p is the chance of gettings "Heads." Half by default

    Returns a RandomSymbol.

    >>> from sympy.stats import Coin, Density
    >>> from sympy import Rational

    >>> X = Coin() # A fair coin toss
    >>> Density(X)
    {H: 1/2, T: 1/2}

    >>> X = Coin(Rational(3, 5)) # An unfair coin
    >>> Density(X)
    {H: 3/5, T: 2/5}
    """
    return CoinPSpace(p, symbol).value

class BinomialPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a binomial distribution.

    This class is for internal use.

    Create Binomial Random Symbols using Binomial function.

    Examples
    ========

    >>> from sympy.stats import Binomial, Density
    >>> from sympy import S

    >>> X = Binomial(4, S.Half) # Four "coin flips"
    >>> Density(X)
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}
    """

    def __new__(cls, n, p, succ=1, fail=0, symbol=None):
        n, p, succ, fail = map(sympify, (n, p, succ, fail))
        density = dict((k*succ + (n-k)*fail,
                binomial(n, k) * p**k * (1-p)**(n-k)) for k in range(0, n+1))
        return create_SingleFinitePSpace(density, symbol, cls)

def Binomial(n, p, succ=1, fail=0, symbol=None):
    """
    Create a Finite Random Variable representing a binomial distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Binomial, Density
    >>> from sympy import S

    >>> X = Binomial(4, S.Half) # Four "coin flips"
    >>> Density(X)
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}
    """

    return BinomialPSpace(n, p, succ, fail, symbol).value

class HypergeometricPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a hypergeometric distribution.

    This class is for internal use.

    Create Hypergeometric Random Symbols using Hypergeometric function.

    Examples
    ========

    >>> from sympy.stats import Hypergeometric, Density
    >>> from sympy import S

    >>> X = Hypergeometric(10, 5, 3) # 10 marbles, 5 white (success), 3 draws
    >>> Density(X)
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}
    """

    def __new__(cls, N, m, n, symbol=None):
        N, m, n = map(sympify, (N, m, n))
        density = dict((k, binomial(m, k) * binomial(N-m, n-k) / binomial(N, n))
                for k in range(max(0, n+m-N), min(m, n) + 1))
        return create_SingleFinitePSpace(density, symbol, cls)

def Hypergeometric(N, m, n, symbol=None):
    """
    Create a Finite Random Variable representing a hypergeometric distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Hypergeometric, Density
    >>> from sympy import S

    >>> X = Hypergeometric(10, 5, 3) # 10 marbles, 5 white (success), 3 draws
    >>> Density(X)
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}
    """

    return HypergeometricPSpace(N, m, n, symbol).value
