"""Sparse rational function fields. """

from sympy.core.sympify import CantSympify
from sympy.polys.monomialtools import lex
from sympy.polys.polyerrors import ExactQuotientFailed
from sympy.polys.rings import PolyRing

def field(sgens, domain, order=lex):
    """Construct new rational function field returning (field, x1, ..., xn). """
    _ring = FracField(sgens, domain, order)
    return (_ring,) + _ring.gens

def xfield(sgens, domain, order=lex):
    """Construct new rational function field returning (field, (x1, ..., xn)). """
    _ring = FracField(sgens, domain, order)
    return (_ring, _ring.gens)

class FracField(object):
    def __init__(self, sgens, domain, order):
        self.ring = PolyRing(sgens, domain, order)
        self.gens = self._gens()

    def _gens(self):
        """Return a list of polynomial generators. """
        return tuple([ FracElement(self, gen) for gen in self.ring.gens ])

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "Rational function field in %s over %s with %s order" % (", ".join(self.ring.sgens), self.ring.domain, self.ring.order)

class FracElement(CantSympify):
    """Sparse rational function. """

    def __init__(self, field, numer, denom=None):
        if denom is not None:
            if denom == numer.ring.zero:
                raise ExactQuotientFailed(numer, denom, numer.ring.domain)
        else:
            denom = numer.ring.one

        self.field = field
        self.numer = field.ring(numer)
        self.denom = field.ring(denom)

    def raw_new(f, numer, denom):
        return FracElement(f.field, numer, denom)
    def new(f, numer, denom):
        return f.raw_new(*numer.cancel(denom))

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        n, d = self.numer, self.denom
        sn, sd = n.str(False), d.str(False)

        if d == d.ring.one:
            return sn

        if len(n) > 1:
            sn = "(" + sn + ")"
        if len(d) > 1:
            sd = "(" + sd + ")"

        return sn + "/" + sd

    def __eq__(f, g):
        return f.numer == g.numer and f.denom == g.denom

    def __bool__(f):
        return bool(f.numer)

    __nonzero__ = __bool__

    def __ne__(f, g):
        return not f.__eq__(g)

    def __pos__(f):
        """Negate all cefficients in ``f``. """
        return f.raw_new(f.numer, f.denom)

    def __neg__(f):
        """Negate all cefficients in ``f``. """
        return f.raw_new(-f.numer, f.denom)

    def __add__(f, g):
        """Add rational functions ``f`` and ``g``. """
        if isinstance(g, FracElement):
            numer = f.numer*g.denom + f.denom*g.numer
            denom = f.denom*g.denom
        else:
            numer = f.numer + f.denom*g
            denom = f.denom

        return f.new(numer, denom)

    def __radd__(f, c):
        numer = f.numer + f.denom*c
        denom = f.denom
        return f.new(numer, denom)

    def __sub__(f, g):
        """Subtract rational functions ``f`` and ``g``. """
        if isinstance(g, FracElement):
            numer = f.numer*g.denom - f.denom*g.numer
            denom = f.denom*g.denom
        else:
            numer = f.numer - f.denom*g
            denom = f.denom

        return f.new(numer, denom)

    def __rsub__(f, c):
        numer = -f.numer + f.denom*c
        denom =  f.denom
        return f.new(numer, denom)

    def __mul__(f, g):
        """Multiply rational functions ``f`` and ``g``. """
        if isinstance(g, FracElement):
            return f.new(f.numer*g.numer, f.denom*g.denom)
        else:
            return f.new(f.numer*g, f.denom)

    def __rmul__(f, c):
        return f.new(f.numer*c, f.denom)

    def __truediv__(f, g):
        """Computes quotient of fractions ``f`` and ``g``. """
        if isinstance(g, FracElement):
            return f.new(f.numer*g.denom, f.denom*g.numer)
        else:
            return f.new(f.numer, f.denom*g)

    __div__ = __truediv__

    def __rtruediv__(f, c):
        return f.new(f.denom*c, f.numer)

    __rdiv__ = __rtruediv__

    def __pow__(f, n):
        """Raise ``f`` to a non-negative power ``n``. """
        if n >= 0:
            return f.raw_new(f.numer**n, f.denom**n)
        else:
            return f.raw_new(f.denom**-n, f.numer**-n)
