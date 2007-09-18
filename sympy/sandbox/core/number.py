
from basic import Basic, Atom
from methods import ArithMeths

class Number(ArithMeths, Atom):

    def __new__(cls, x, **options):
        if isinstance(x, Basic): return x
        if isinstance(x, (int, long)): return Integer(x, **options)
        return Basic.__new__(cls, x, **options)
        
class Real(Number):

    pass

class Rational(Number):

    def __new__(cls, p, q=1):
        return Fraction(p, q)

class Fraction(Rational, tuple):

    def __new__(cls, p, q, **options):
        assert not options,`options`
        if q==1:
            return Integer(p)
        assert isinstance(p, int),`p`
        assert isinstance(q, int),`q`
        return tuple.__new__(cls, (p,q))

    @property
    def p(self): return self[0]
    @property
    def q(self): return self[1]

    def torepr(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.p, self.q)

class Integer(Rational, int):

    def __new__(cls, p, **options):
        assert not options,`options`
        obj = int.__new__(cls, p)
        return obj

    @property
    def p(self): return int(self)
    @property
    def q(self): return 1

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.p)

    def __add__(self, other):
        if isinstance(other, int):
            return Integer(int(self) + int(other))
        return Basic.Add(self, other)
