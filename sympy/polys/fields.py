"""Sparse rational function fields. """

from operator import add, mul

from sympy.core.expr import Expr
from sympy.core.sympify import CantSympify, sympify
from sympy.polys.rings import PolyElement
from sympy.polys.monomialtools import lex
from sympy.polys.polyerrors import ExactQuotientFailed, CoercionFailed
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.fractionfield import FractionField
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
        try:
            return self.new(self.ring.ground_new(element))
        except CoercionFailed:
            domain = self.domain

            if domain.has_Ring and domain.has_assoc_Field:
                ground_field = domain.get_field()
                element = ground_field.convert(element)
                numer = ground_field.numer(element)
                denom = ground_field.denom(element)
                return self.new(numer, denom)
            else:
                raise

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
        elif isinstance(element, tuple) and len(element) == 2:
            numer, denom = map(self.ring.ring_new, element)
            return self.new(numer, denom)
        elif isinstance(element, basestring):
            raise NotImplementedError("parsing")
        elif isinstance(element, Expr):
            return self.from_expr(element)
        else:
            return self.ground_new(element)

    __call__ = field_new

    def _rebuild_expr(self, expr, mapping):
        domain = self.domain

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
                try:
                    return domain.convert(expr)
                except CoercionFailed:
                    if domain.has_Ring and domain.has_assoc_Field:
                        return domain.get_field().convert(expr)
                    else:
                        raise

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        mapping = dict(zip(self.symbols, self.gens))

        try:
            frac = self._rebuild_expr(expr, mapping)
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
        from sympy.polys.domains.fractionfield import FractionField
        return FractionField(self)

    def to_ring(self):
        from sympy.polys.rings import PolyRing
        return PolyRing(self.symbols, self.domain, self.order)

class FracElement(DomainElement, DefaultPrinting, CantSympify):
    """Sparse rational function. """

    def __init__(self, field, numer, denom=None):
        if denom is not None:
            if not denom:
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

    def parent(self):
        return self.field.to_domain()

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.field, self.numer, self.denom))
        return _hash

    def copy(self):
        return self.raw_new(self.numer.copy(), self.denom.copy())

    def set_field(self, new_field):
        if self.field == new_field:
            return self
        else:
            new_ring = new_field.ring
            numer = self.numer.set_ring(new_ring)
            denom = self.denom.set_ring(new_ring)
            return new_field.new(numer, denom)

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

    def _extract_ground(self, element):
        domain = self.field.domain

        try:
            element = domain.convert(element)
        except CoercionFailed:
            if domain.has_Ring and domain.has_assoc_Field:
                ground_field = domain.get_field()

                try:
                    element = ground_field.convert(element)
                except CoercionFailed:
                    pass
                else:
                    return -1, ground_field.numer(element), ground_field.denom(element)

            return 0, None, None
        else:
            return 1, element, None

    def __add__(f, g):
        """Add rational functions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if f.field == g.field:
                return f.new(f.numer*g.denom + f.denom*g.numer, f.denom*g.denom)
            elif isinstance(field.domain, FractionField) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionField) and g.field.domain.field == field:
                return g.__radd__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement):
            if field.ring == g.ring:
                return f.new(f.numer + f.denom*g, f.denom)
            else:
                return g.__radd__(f)

        return f.__radd__(g)

    def __radd__(f, c):
        if isinstance(c, PolyElement) and f.field.ring == c.ring:
            return f.new(f.numer + f.denom*c, f.denom)

        op, g_numer, g_denom = f._extract_ground(c)

        if op == 1:
            return f.new(f.numer + f.denom*g_numer, f.denom)
        elif not op:
            return NotImplemented
        else:
            return f.new(f.numer*g_denom + f.denom*g_numer, f.denom*g_denom)

    def __sub__(f, g):
        """Subtract rational functions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if f.field == g.field:
                return f.new(f.numer*g.denom - f.denom*g.numer, f.denom*g.denom)
            elif isinstance(field.domain, FractionField) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionField) and g.field.domain.field == field:
                return g.__rsub__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement):
            if field.ring == g.ring:
                return f.new(f.numer - f.denom*g, f.denom)
            else:
                return g.__rsub__(f)

        op, g_numer, g_denom = f._extract_ground(g)

        if op == 1:
            return f.new(f.numer - f.denom*g_numer, f.denom)
        elif not op:
            return NotImplemented
        else:
            return f.new(f.numer*g_denom - f.denom*g_numer, f.denom*g_denom)

    def __rsub__(f, c):
        if isinstance(c, PolyElement) and f.field.ring == c.ring:
            return f.new(-f.numer + f.denom*c, f.denom)

        op, g_numer, g_denom = f._extract_ground(c)

        if op == 1:
            return f.new(-f.numer + f.denom*g_numer, f.denom)
        elif not op:
            return NotImplemented
        else:
            return f.new(-f.numer*g_denom + f.denom*g_numer, f.denom*g_denom)

    def __mul__(f, g):
        """Multiply rational functions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if field == g.field:
                return f.new(f.numer*g.numer, f.denom*g.denom)
            elif isinstance(field.domain, FractionField) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionField) and g.field.domain.field == field:
                return g.__rmul__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement):
            if field.ring == g.ring:
                return f.new(f.numer*g, f.denom)
            else:
                return g.__rmul__(f)

        return f.__rmul__(g)

    def __rmul__(f, c):
        if isinstance(c, PolyElement) and f.field.ring == c.ring:
            return f.new(f.numer*c, f.denom)

        op, g_numer, g_denom = f._extract_ground(c)

        if op == 1:
            return f.new(f.numer*g_numer, f.denom)
        elif not op:
            return NotImplemented
        else:
            return f.new(f.numer*g_numer, f.denom*g_denom)

    def __truediv__(f, g):
        """Computes quotient of fractions ``f`` and ``g``. """
        field = f.field

        if isinstance(g, FracElement):
            if field == g.field:
                return f.new(f.numer*g.denom, f.denom*g.numer)
            elif isinstance(field.domain, FractionField) and field.domain.field == g.field:
                pass
            elif isinstance(g.field.domain, FractionField) and g.field.domain.field == field:
                return g.__rtruediv__(f)
            else:
                return NotImplemented
        elif isinstance(g, PolyElement):
            if field.ring == g.ring:
                return f.new(f.numer, f.denom*g)
            else:
                return g.__rtruediv__(f)

        op, g_numer, g_denom = f._extract_ground(g)

        if op == 1:
            return f.new(f.numer, f.denom*g_numer)
        elif not op:
            return NotImplemented
        else:
            return f.new(f.numer*g_denom, f.denom*g_numer)

    __div__ = __truediv__

    def __rtruediv__(f, c):
        if isinstance(c, PolyElement) and f.field.ring == c.ring:
            return f.new(f.denom*c, f.numer)

        op, g_numer, g_denom = f._extract_ground(c)

        if op == 1:
            return f.new(f.denom*g_numer, f.numer)
        elif not op:
            return NotImplemented
        else:
            return f.new(f.denom*g_numer, f.numer*g_denom)

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
