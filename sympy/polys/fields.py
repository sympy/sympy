"""Sparse rational function fields. """

from sympy.core.sympify import CantSympify
from sympy.polys.monomialtools import lex
from sympy.polys.polyerrors import ExactQuotientFailed
from sympy.polys.rings import PolyRing

def field(symbols, domain, order=lex):
    """Construct new rational function field returning (field, x1, ..., xn). """
    _field = FracField(symbols, domain, order)
    return (_field,) + _field.gens

def xfield(symbols, domain, order=lex):
    """Construct new rational function field returning (field, (x1, ..., xn)). """
    _field = FracField(symbols, domain, order)
    return (_field, _field.gens)

def vfield(symbols, domain, order=lex):
    """Construct new rational function field and inject generators into global namespace. """
    from inspect import currentframe
    frame = currentframe().f_back

    try:
        _field = FracField(symbols, domain, order)

        for name, gen in zip(_field.symbols, _field.gens):
            frame.f_globals[name] = gen
    finally:
        del frame  # break cyclic dependencies as stated in inspect docs

    return (_field, _field.gens)

class FracField(object):
    def __init__(self, symbols, domain, order):
        self.ring = PolyRing(symbols, domain, order)
        self.gens = self._gens()

    def _gens(self):
        """Return a list of polynomial generators. """
        return tuple([ FracElement(self, gen) for gen in self.ring.gens ])

    def __repr__(self):
        return "%s(%s, %s, %s)" % (self.__class__.__name__, repr(self.ring.symbols), repr(self.ring.domain), repr(self.ring.order))

    def __str__(self):
        return "Rational function field in %s over %s with %s order" % (", ".join(map(str, self.ring.symbols)), self.ring.domain, self.ring.order)

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

    def to_poly(f):
        assert f.denom == 1
        return f.numer

    def __repr__(self):
        numer_terms = list(self.numer.terms())
        numer_terms.sort(key=self.field.ring.order, reverse=True)
        denom_terms = list(self.denom.terms())
        denom_terms.sort(key=self.field.ring.order, reverse=True)
        return "%s(%s, %s, %s)" % (self.__class__.__name__, repr(self.field), repr(numer_terms), repr(denom_terms))

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

    def diff(f, x):
        if isinstance(x, list) and a is None:
            x = [ X.to_poly() for X in x ]
        else:
            x = x.to_poly()
        return f.new(f.numer.diff(x)*f.denom - f.numer*f.denom.diff(x), f.denom**2)

    def subs(f, x, a=None):
        if isinstance(x, list) and a is None:
            x = [ (X.to_poly(), a) for X, a in x ]
            return f.new(f.numer.subs(x), f.denom.subs(x))
        else:
            x = x.to_poly()
            return f.new(f.numer.subs(x, a), f.denom.subs(x, a))
