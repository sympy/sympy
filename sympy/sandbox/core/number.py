
from basic import Atom
from methods import ArithMeths

class Number(ArithMeths, Atom):
    pass

class Real(Number):

    pass

class Rational(Number):

    def __new__(cls, p, q=1, **options):
        assert not options,`options`
        if q==1:
            return Integer(p)
        assert isinstance(p, int),`p`
        assert isinstance(q, int),`q`
        obj = object.__new__(cls)
        obj.p = int(p)
        obj.q = int(q)
        return obj

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

    def __new__(cls, p, **options):
        assert not options,`options`
        obj = int.__new__(cls, p)
        obj.p = int(obj)
        obj.q = 1
        return obj

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.p)
