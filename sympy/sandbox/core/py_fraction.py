from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational

class Fraction(Rational, tuple):

    @memoizer_immutable_args('Fraction.__new__')
    def __new__(cls, p, q):
        p, q = sympify(p), sympify(q)
        if q==1: return p
        if p.is_Integer and q.is_Integer:
            if q.is_negative:
                p,q = -p,-q
            return tuple.__new__(cls, (p,q))
        if q.is_Fraction:
            iq = Fraction(q.q, q.p)
            return p * iq
        if p.is_Fraction:
            return p * q
        raise TypeError(`p,q`)

    @property
    def p(self): return self[0]
    
    @property
    def q(self): return self[1]

    def __int__(self):
        raise NotImplementedError
