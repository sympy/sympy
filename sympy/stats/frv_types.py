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

from sympy.stats.frv import SingleFinitePSpace
from sympy import S, sympify, Rational, binomial

__all__ = ['FiniteRV', 'DiscreteUniform', 'Die', 'Bernoulli', 'Coin',
        'Binomial', 'Hypergeometric']

def FiniteRV(name, density):
    """
    Create a Finite Random Variable given a dict representing the density.

    Returns a RandomSymbol.

    >>> from sympy.stats import FiniteRV, P, E

    >>> density = {0: .1, 1: .2, 2: .3, 3: .4}
    >>> X = FiniteRV('X', density)

    >>> E(X)
    2.00000000000000
    >>> P(X>=2)
    0.700000000000000
    """
    return SingleFinitePSpace.fromdict(name, density).value

class DiscreteUniformPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a discrete uniform
    distribution.

    This class is for internal use.

    Create DiscreteUniform Random Symbols using DiscreteUniform function

    Examples
    ========

    >>> from sympy.stats import DiscreteUniform, density
    >>> from sympy import symbols

    >>> X = DiscreteUniform('X', symbols('a b c')) # equally likely over a, b, c
    >>> density(X)
    {a: 1/3, b: 1/3, c: 1/3}

    >>> Y = DiscreteUniform('Y', range(5)) # distribution over a range
    >>> density(Y)
    {0: 1/5, 1: 1/5, 2: 1/5, 3: 1/5, 4: 1/5}
    """
    def __new__(cls, name, items):
        density = dict((sympify(item), Rational(1, len(items)))
                       for item in items)
        return cls.fromdict(name, density)

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
    >>> density(X)
    {a: 1/3, b: 1/3, c: 1/3}

    >>> Y = DiscreteUniform('Y', range(5)) # distribution over a range
    >>> density(Y)
    {0: 1/5, 1: 1/5, 2: 1/5, 3: 1/5, 4: 1/5}

    """
    return DiscreteUniformPSpace(name, items).value

class DiePSpace(DiscreteUniformPSpace):
    """
    Create a Finite Random Variable representing a fair die.

    This class is for internal use.

    Create Dice Random Symbols using Die function

    >>> from sympy.stats import Die, density

    >>> D6 = Die('D6', 6) # Six sided Die
    >>> density(D6)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}

    >>> D4 = Die('D4', 4) # Four sided Die
    >>> density(D4)
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}
    """
    def __new__(cls, name, sides):
        return DiscreteUniformPSpace.__new__(cls, name, range(1, sides+1))

def Die(name, sides=6):
    """
    Create a Finite Random Variable representing a fair die.

    Returns a RandomSymbol.

    >>> from sympy.stats import Die, density

    >>> D6 = Die('D6', 6) # Six sided Die
    >>> density(D6)
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}

    >>> D4 = Die('D4', 4) # Four sided Die
    >>> density(D4)
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}
    """

    return DiePSpace(name, sides).value

class BernoulliPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a Bernoulli process.

    Returns a RandomSymbol.

    This class is for internal use.

    Create Bernoulli Random Symbols using Bernoulli function.

    >>> from sympy.stats import Bernoulli, density
    >>> from sympy import S

    >>> X = Bernoulli('X', S(3)/4) # 1-0 Bernoulli variable, probability = 3/4
    >>> density(X)
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli('X', S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> density(X)
    {Heads: 1/2, Tails: 1/2}
    """

    def __new__(cls, name,  p, succ, fail):
        succ, fail, p = map(sympify, (succ, fail, p))
        density = {succ: p, fail: (1-p)}
        return cls.fromdict(name, density)

def Bernoulli(name, p, succ=1, fail=0):
    """
    Create a Finite Random Variable representing a Bernoulli process.

    Returns a RandomSymbol

    >>> from sympy.stats import Bernoulli, density
    >>> from sympy import S

    >>> X = Bernoulli('X', S(3)/4) # 1-0 Bernoulli variable, probability = 3/4
    >>> density(X)
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli('X', S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> density(X)
    {Heads: 1/2, Tails: 1/2}
    """

    return BernoulliPSpace(name, p, succ, fail).value

class CoinPSpace(BernoulliPSpace):
    """
    A probability space representing a coin toss.

    Probability p is the chance of gettings "Heads." Half by default

    This class is for internal use.

    Create Coin's using Coin function

    >>> from sympy.stats import Coin, density
    >>> from sympy import Rational

    >>> C = Coin('C') # A fair coin toss
    >>> density(C)
    {H: 1/2, T: 1/2}

    >>> C2 = Coin('C2', Rational(3, 5)) # An unfair coin
    >>> density(C2)
    {H: 3/5, T: 2/5}
    """
    def __new__(cls, name, p):
        return BernoulliPSpace.__new__(cls, name, p, 'H', 'T')

def Coin(name, p=S.Half):
    """
    Create a Finite Random Variable representing a Coin toss.

    Probability p is the chance of gettings "Heads." Half by default

    Returns a RandomSymbol.

    >>> from sympy.stats import Coin, density
    >>> from sympy import Rational

    >>> C = Coin('C') # A fair coin toss
    >>> density(C)
    {H: 1/2, T: 1/2}

    >>> C2 = Coin('C2', Rational(3, 5)) # An unfair coin
    >>> density(C2)
    {H: 3/5, T: 2/5}
    """
    return CoinPSpace(name, p).value

class BinomialPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a binomial distribution.

    This class is for internal use.

    Create Binomial Random Symbols using Binomial function.

    Examples
    ========

    >>> from sympy.stats import Binomial, density
    >>> from sympy import S

    >>> X = Binomial('X', 4, S.Half) # Four "coin flips"
    >>> density(X)
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}
    """

    def __new__(cls, name, n, p, succ, fail):
        n, p, succ, fail = map(sympify, (n, p, succ, fail))
        density = dict((k*succ + (n-k)*fail,
                binomial(n, k) * p**k * (1-p)**(n-k)) for k in range(0, n+1))
        return cls.fromdict(name, density)

def Binomial(name, n, p, succ=1, fail=0):
    """
    Create a Finite Random Variable representing a binomial distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Binomial, density
    >>> from sympy import S

    >>> X = Binomial('X', 4, S.Half) # Four "coin flips"
    >>> density(X)
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}
    """

    return BinomialPSpace(name, n, p, succ, fail).value

class HypergeometricPSpace(SingleFinitePSpace):
    """
    Create a Finite Random Variable representing a hypergeometric distribution.

    This class is for internal use.

    Create Hypergeometric Random Symbols using Hypergeometric function.

    Examples
    ========

    >>> from sympy.stats import Hypergeometric, density
    >>> from sympy import S

    >>> X = Hypergeometric('X', 10, 5, 3) # 10 marbles, 5 white (success), 3 draws
    >>> density(X)
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}
    """

    def __new__(cls, name, N, m, n):
        N, m, n = map(sympify, (N, m, n))
        density = dict((k, binomial(m, k) * binomial(N-m, n-k) / binomial(N, n))
                for k in range(max(0, n+m-N), min(m, n) + 1))
        return cls.fromdict(name, density)

def Hypergeometric(name, N, m, n):
    """
    Create a Finite Random Variable representing a hypergeometric distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Hypergeometric, density
    >>> from sympy import S

    >>> X = Hypergeometric('X', 10, 5, 3) # 10 marbles, 5 white (success), 3 draws
    >>> density(X)
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}
    """

    return HypergeometricPSpace(name, N, m, n).value
