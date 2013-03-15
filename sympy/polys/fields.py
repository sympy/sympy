"""Sparse rational function fields. """

from sympy.core.expr import Expr
from sympy.core.sympify import CantSympify, sympify
from sympy.polys.rings import PolyElement
from sympy.polys.monomialtools import lex
from sympy.polys.polyerrors import ExactQuotientFailed, CoercionFailed
from sympy.polys.domains.fractionfield import FractionFieldNG
from sympy.printing.defaults import DefaultPrinting

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

        for sym, gen in zip(_field.symbols, _field.gens):
            frame.f_globals[sym.name] = gen
    finally:
        del frame  # break cyclic dependencies as stated in inspect docs

    return _field

class FracField(DefaultPrinting):

    def __init__(self, symbols, domain, order):
        from sympy.polys.rings import PolyRing
        self.ring = PolyRing(symbols, domain, order)
        self.dtype = FracElement
        self.symbols = self.ring.symbols
        self.ngens = len(self.symbols)
        self.domain = self.ring.domain
        self.order = self.ring.order
        self.gens = self._gens()

    def _gens(self):
        """Return a list of polynomial generators. """
        return tuple([ self.dtype(self, gen) for gen in self.ring.gens ])

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.symbols, self.domain, self.order))
        return _hash

    def __eq__(self, other):
        return isinstance(other, FracField) and self.ring == other.ring

    def __ne__(self, other):
        return not self.__eq__(other)

    def new(self, numer, denom=None):
        return self.dtype(self, numer, denom)

    def domain_new(self, element):
        return self.domain.convert(element)

    def ground_new(self, element):
        return self.new(self.ring.ground_new(element))

    def field_new(self, element):
        if isinstance(element, FracElement):
            if self == element.field:
                return element
            else:
                raise NotImplementedError("conversion")
        elif isinstance(element, PolyElement):
            if self.ring == element.ring:
                return self.new(element)
            else:
                raise NotImplementedError("conversion")
        elif isinstance(element, basestring):
            raise NotImplementedError("parsing")
        elif isinstance(element, Expr):
            return self.from_expr(element)
        else:
            return self.ground_new(element)

    __call__ = field_new

    @staticmethod
    def _rebuild_expr(expr, mapping, domain):
        from operator import add, mul

        def _rebuild(expr):
            generator = mapping.get(expr)

            if generator is not None:
                return generator
            elif expr.is_Add:
                return reduce(add, map(_rebuild, expr.args))
            elif expr.is_Mul:
                return reduce(mul, map(_rebuild, expr.args))
            elif expr.is_Pow and expr.exp.is_Integer:
                return _rebuild(expr.base)**int(expr.exp)
            else:
                return domain.convert(expr)

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        mapping = dict(zip(self.symbols, self.gens))

        try:
            frac = self._rebuild_expr(expr, mapping, self.domain)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a rational function in %s, got %s" % (self, expr))
        else:
            return self.field_new(frac)

    @property
    def zero(self):
        return self.new(self.ring.zero)

    @property
    def one(self):
        return self.new(self.ring.one)

    def to_domain(self):
        from sympy.polys.domains.fractionfield import FractionFieldNG
        return FractionFieldNG(self)

    def to_ring(self):
        from sympy.polys.rings import PolyRing
        return PolyRing(self.symbols, self.domain, self.order)

class FracElement(CantSympify, DefaultPrinting):
    """Sparse rational function. """

    def __init__(self, field, numer, denom=None):
        if denom is not None:
            if denom == numer.ring.zero:
                raise ZeroDivisionError
        else:
            denom = numer.ring.one

        self.field = field
        self.numer = field.ring(numer)
        self.denom = field.ring(denom)

    def raw_new(f, numer, denom):
        return f.__class__(f.field, numer, denom)
    def new(f, numer, denom):
        return f.raw_new(*numer.cancel(denom))

    def to_poly(f):
        assert f.denom == 1
        return f.numer

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.field, self.numer, self.denom))
        return _hash

    def copy(self):
        return self.raw_new(self.numer.copy(), self.denom.copy())

    def as_expr(self, *symbols):
        return self.numer.as_expr(*symbols)/self.denom.as_expr(*symbols)

    def __eq__(f, g):
        if isinstance(g, FracElement):
            return f.numer == g.numer and f.denom == g.denom
        else:
            return f.numer == g and f.denom == 1

    def __ne__(f, g):
        return not f.__eq__(g)

    def __bool__(f):
        return bool(f.numer)

    __nonzero__ = __bool__

    def __pos__(f):
        """Negate all cefficients in ``f``. """
        return f.raw_new(f.numer, f.denom)

    def __neg__(f):
        """Negate all cefficients in ``f``. """
        return f.raw_new(-f.numer, f.denom)

    def __add__(f, g):
        """Add rational functions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if f.field == g.field:
                return f.new(f.numer*g.denom + f.denom*g.numer, f.denom*g.denom)
            elif isinstance(field.domain, FractionFieldNG) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionFieldNG) and g.field.domain.field == field:
                return g.__radd__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement) and field.ring != g.ring:
            return g.__radd__(f)

        return f.new(f.numer + f.denom*g, f.denom)

    def __radd__(f, c):
        numer = f.numer + f.denom*c
        denom = f.denom
        return f.new(numer, denom)

    def __sub__(f, g):
        """Subtract rational functions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if f.field == g.field:
                return f.new(f.numer*g.denom - f.denom*g.numer, f.denom*g.denom)
            elif isinstance(field.domain, FractionFieldNG) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionFieldNG) and g.field.domain.field == field:
                return g.__rsub__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement) and field.ring != g.ring:
            return g.__rsub__(f)

        return f.new(f.numer - f.denom*g, f.denom)

    def __rsub__(f, c):
        numer = -f.numer + f.denom*c
        denom =  f.denom
        return f.new(numer, denom)

    def __mul__(f, g):
        """Multiply rational functions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if field == g.field:
                return f.new(f.numer*g.numer, f.denom*g.denom)
            elif isinstance(field.domain, FractionFieldNG) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionFieldNG) and g.field.domain.field == field:
                return g.__rmul__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement) and field.ring != g.ring:
            return g.__rmul__(f)

        return f.new(f.numer*g, f.denom)

    def __rmul__(f, c):
        return f.new(f.numer*c, f.denom)

    def __truediv__(f, g):
        """Computes quotient of fractions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if field == g.field:
                return f.new(f.numer*g.denom, f.denom*g.numer)
            elif isinstance(field.domain, FractionFieldNG) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionFieldNG) and g.field.domain.field == field:
                return g.__rtruediv__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement) and field.ring != g.ring:
            return g.__rtruediv__(f)

        return f.new(f.numer, f.denom*g)

    __div__ = __truediv__

    def __rtruediv__(f, c):
        return f.new(f.denom*c, f.numer)

    __rdiv__ = __rtruediv__

    def __pow__(f, n):
        """Raise ``f`` to a non-negative power ``n``. """
        if isinstance(n, int):
            if n >= 0:
                return f.raw_new(f.numer**n, f.denom**n)
            else:
                return f.raw_new(f.denom**-n, f.numer**-n)
        else:
            return NotImplemented

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
