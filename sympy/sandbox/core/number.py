
from basic import Atom
from methods import ArithMeths

class Number(ArithMeths, Atom):
    pass

class Real(Number):

    pass

class Rational(Number):

    @classmethod
    def canonize(cls, p, q=1):
        if q==1:
            return Integer(p)
        return (), dict(p=p,q=q)

    def __hash__(self):
        try:
            return self._hash
        except AttributeError:
            h = hash((self.__class__.__name__, self.p, self.q))
            self._hash = h
            return h

    def torepr(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.p, self.q)

class Integer(Rational, int):

    _new = int.__new__
    
    @classmethod
    def canonize(cls, p):
        p = int(p)
        return (p,), {}

    @property
    def p(self): return int(self)

    @property
    def q(self): return 1

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.p)
