from sympy.statistics.frv import SingleFinitePSpace, create_SingleFinitePSpace
from sympy import S, sympify, Rational

def FiniteRV(density, symbol=None):
    """
    Create a Finite Random Variable given a dict representing the density
    Returns a RandomSymbol

    >>> from sympy.statistics import FiniteRV, P, E

    >>> density = {0: .1, 1: .2, 2: .3, 3: .4}
    >>> X = FiniteRV(density)

    >>> E(X)
    2.00000000000000
    >>> P(X>=2)
    0.700000000000000


    """
    return create_SingleFinitePSpace(density, symbol).value

class DiePSpace(SingleFinitePSpace):
    """
    Create a Finite Random Varible representing a fair die

    This class is for internal use.

    Create Dice Random Symbols using Die function

    >>> from sympy.statistics import Die, Density

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
        density = dict((i,Rational(1,sides)) for i in range(1,sides+1))
        return create_SingleFinitePSpace(density, symbol, cls)

def Die(sides=6, symbol=None):
    """
    Create a Finite Random Varible representing a fair die
    Returns a RandomSymbol

    >>> from sympy.statistics import Die, Density

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
    Create a Finite Random Varible representing a Bernoulli process
    Returns a RandomSymbol

    This class is for internal use.

    Create Bernoulli Random Symbols using Bernoulli function

    >>> from sympy.statistics import Bernoulli, Density
    >>> from sympy import S

    >>> X = Bernoulli(S(3)/4, 1, 0) # 1-0 Bernoulli variable, probability = 3/4
    >>> Density(X)
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli(S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> Density(X)
    {Heads: 1/2, Tails: 1/2}

    """

    _count = 0
    _name = 'bernoulli'
    def __new__(cls, p, a, b, symbol=None):
        a, b, p = map(sympify, (a, b, p))
        density = {a:p, b:(1-p)}
        return create_SingleFinitePSpace(density, symbol, cls)

def Bernoulli(p, a, b, symbol=None):
    """
    Create a Finite Random Varible representing a Bernoulli process
    Returns a RandomSymbol

    >>> from sympy.statistics import Bernoulli, Density
    >>> from sympy import S

    >>> X = Bernoulli(S(3)/4, 1, 0) # 1-0 Bernoulli variable, probability = 3/4
    >>> Density(X)
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli(S.Half, 'Heads', 'Tails') # A fair coin toss
    >>> Density(X)
    {Heads: 1/2, Tails: 1/2}

    """

    return BernoulliPSpace(p, a, b, symbol).value

class CoinPSpace(BernoulliPSpace):
    """
    A probability space representing a coin toss
    Probability p is the chance of gettings "Heads." Half by default

    This class is for internal use.

    Create Coin's using Coin function

    >>> from sympy.statistics import Coin, Density
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
    Create a Finite Random Varible representing a Coin toss
    Probability p is the chance of gettings "Heads." Half by default

    Returns a RandomSymbol

    >>> from sympy.statistics import Coin, Density
    >>> from sympy import Rational

    >>> X = Coin() # A fair coin toss
    >>> Density(X)
    {H: 1/2, T: 1/2}

    >>> X = Coin(Rational(3, 5)) # An unfair coin
    >>> Density(X)
    {H: 3/5, T: 2/5}

    """
    return CoinPSpace(p, symbol).value



