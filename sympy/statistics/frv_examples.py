from sympy.statistics.frv import SingleFinitePSpace, create_SingleFinitePSpace
from sympy import S, sympify, Rational

class DiePSpace(SingleFinitePSpace):
    _count = 0
    _name = 'die'
    def __new__(cls, sides=6, symbol=None):
        density = dict((i,Rational(1,sides)) for i in range(1,sides+1))
        return create_SingleFinitePSpace(density, symbol, cls)

def Die(sides=6, symbol=None):
    return DiePSpace(sides, symbol).value

class BernoulliPSpace(SingleFinitePSpace):
    _count = 0
    _name = 'bernoulli'
    def __new__(cls, p, a, b, symbol=None):
        a, b, p = map(sympify, (a, b, p))
        density = {a:p, b:(1-p)}
        return create_SingleFinitePSpace(density, symbol, cls)

def Bernoulli(p, a, b, symbol=None):
    return BernoulliPSpace(p, a, b, symbol).value

class CoinPSpace(BernoulliPSpace):
    _count = 0
    _name = 'coin'
    def __new__(cls, p=S.Half, symbol=None):
        return Bernoulli.__new__(cls, p, 'H', 'T', symbol)

def Coin(p=S.Half, symbol=None):
    return CoinPSpace(p, symbol).value



