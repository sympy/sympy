from __future__ import annotations

from functools import reduce
from operator import add, mul

from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify, CantSympify
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.laurentpolynomialring import LaurentPolynomialRing
from sympy.polys.fields import FracField
from sympy.polys.orderings import LexOrder
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.rings import PolyElement, PolyRing
from sympy.printing.defaults import DefaultPrinting

from typing import Any, Optional


def laurent_ring(symbols, domain, order='lex'):
    """Construct a Laurent polynomial ring.

    Examples
    ========

    >>> from sympy.polys.laurent import laurent_ring
    >>> from sympy import ZZ
    >>> R, x, y = laurent_ring('x,y', ZZ, order='lex')
    >>> R
    Laurent polynomial ring in x, y over ZZ with lex order
    >>> p = x**2 + 1/y
    >>> p
    x**2 + 1/y
    >>> p ** 2
    x**4 + 2*x**2/y + 1/y**2
    """
    ring = LaurentPolyRing(symbols, domain, order=order)
    return (ring,) + ring.gens


class LaurentPolyRing:
    """Ring of Laurent polynomials.

    Examples
    ========

    >>> from sympy.polys.laurent import LaurentPolyRing
    >>> from sympy import ZZ
    >>> from sympy.abc import x, y
    >>> R = LaurentPolyRing([x, y], ZZ, order='lex')
    >>> R
    Laurent polynomial ring in x, y over ZZ with lex order
    >>> x, y = R.gens
    >>> p = 1/x + x*y
    >>> p
    x*y + 1/x
    >>> p ** 2
    x**2*y**2 + 2*y + 1/x**2
    """
    def __init__(self, symbols, domain, order=LexOrder()):
        numer_ring = PolyRing(symbols, domain, order=order)

        self.numer_ring = numer_ring
        self.domain = numer_ring.domain
        self.symbols = numer_ring.symbols
        self.order = numer_ring.order

        self.zero = self.from_polyelement(numer_ring.zero)
        self.one = self.from_polyelement(numer_ring.one)
        self.gens = tuple(self.from_polyelement(g) for g in numer_ring.gens)

        # Make the symbols accessible as R.x, R.y etc.
        for symbol, generator in zip(self.symbols, self.gens):
            if isinstance(symbol, Symbol):
                name = symbol.name
                if not hasattr(self, name):
                    setattr(self, name, generator)

    def __eq__(self, other):
        if not isinstance(other, LaurentPolyRing):
            return NotImplemented
        return self.numer_ring == other.numer_ring

    def __hash__(self):
        return hash(self.numer_ring)

    def __repr__(self):
        syms = ', '.join(map(str, self.symbols))
        return f"Laurent polynomial ring in {syms} over {self.domain} with {self.order} order"

    def __call__(self, element):
        return self.ring_new(element)

    def new(self, numer, denom=None):
        """Construct a new Laurent polynomial from a numerator and denominator."""
        if denom is None:
            denom = self.numer_ring.one
        return LaurentPolyElement(self, numer, denom)

    def ring_new(self, element):
        if isinstance(element, tuple) and len(element) == 2:
            numer, denom = element
            numer = self.numer_ring.ring_new(numer)
            denom = self.numer_ring.ring_new(denom)
            return self.new(numer, denom)
        elif isinstance(element, Expr):
            return self.from_expr(element)
        else:
            numer = self.numer_ring.ring_new(element)
            return self.new(numer)

    def ground_new(self, coeff):
        """Construct a new Laurent polynomial from an element of the ground domain."""
        return self.from_polyelement(self.numer_ring.ground_new(coeff))

    def from_polyelement(self, numer: PolyElement):
        """Construct a new Laurent polynomial from a polynomial."""
        assert self.numer_ring == numer.ring
        return LaurentPolyElement(self, numer, self.numer_ring.one)

    def from_expr(self, expr):
        """Construct a new Laurent polynomial from a SymPy expression."""
        mapping = dict(list(zip(self.symbols, self.gens)))

        try:
            poly = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a polynomial in %s, got %s" % (self, expr))
        else:
            return poly

    def _rebuild_expr(self, expr, mapping):
        # XXX: This is mostly just a copy from PolyElement. Maybe it could be
        #      refactored into something reusable...
        domain = self.domain

        def _rebuild(expr):
            generator = mapping.get(expr)

            if generator is not None:
                return generator
            elif expr.is_Add:
                return reduce(add, list(map(_rebuild, expr.args)))
            elif expr.is_Mul:
                return reduce(mul, list(map(_rebuild, expr.args)))
            else:
                # XXX: Use as_base_exp() to handle Pow(x, n) and also exp(n)
                # XXX: E can be a generator e.g. sring([exp(2)]) -> ZZ[E]
                base, exp = expr.as_base_exp()
                if exp.is_Integer and exp != 1:
                    return _rebuild(base)**int(exp)
                else:
                    return self.ground_new(domain.convert(expr))

        return _rebuild(sympify(expr))

    def from_dict(self, element, orig_domain=None):
        """Construct a new Laurent polynomial from a dictionary."""
        min_degrees = tuple(map(min, zip(*element.keys())))

        if all(d >= 0 for d in min_degrees):
            numer = self.numer_ring.from_dict(element, orig_domain)
            return self.new(numer)

        # Factor out the denominator monomial
        denom_monom = tuple(max(0, -d) for d in min_degrees)
        monomial_mul = self.numer_ring.monomial_mul
        element = {monomial_mul(m, denom_monom): c for m, c in element.items()}
        numer = self.numer_ring.from_dict(element, orig_domain)
        denom = self.numer_ring.term_new(denom_monom, self.numer_ring.domain.one)

        return self.new(numer, denom)

    def to_field(self) -> FracField:
        """Return a ``FracField`` with the same generators."""
        return FracField(self.symbols, self.domain, self.order)

    def to_domain(self) -> LaurentPolynomialRing:
        return LaurentPolynomialRing(self)


class LaurentPolyElement(DomainElement, DefaultPrinting, CantSympify):
    """A class for representing Laurent polynomials.

    Examples
    ========

    >>> from sympy.polys.laurent import laurent_ring
    >>> from sympy import QQ
    >>> R, x, y = laurent_ring('x,y', QQ)
    >>> p = 1/x + x*y
    >>> p
    x*y + 1/x
    >>> type(p)
    <class 'sympy.polys.laurent.LaurentPolyElement'>
    """
    def __init__(self,
                 ring: LaurentPolyRing,
                 numer: PolyElement,
                 denom: PolyElement
                 ):
        self._check(ring, numer, denom)
        self.ring = ring
        self.numer = numer
        self.denom = denom

    def parent(self) -> LaurentPolynomialRing:
        return self.ring.to_domain()

    def __eq__(self, other):
        if not isinstance(other, LaurentPolyElement):
            return NotImplemented
        return (self.ring == other.ring
                and self.numer == other.numer
                and self.denom == other.denom)

    def __hash__(self):
        return hash((self.ring, self.numer, self.denom))

    @classmethod
    def _check(self,
               ring: LaurentPolyRing,
               numer: PolyElement,
               denom: PolyElement,
               check_cancelled: Optional[bool] = True
               ) -> None:
        """Validate the internal invariants."""
        assert type(ring) is LaurentPolyRing
        assert isinstance(numer, PolyElement)
        assert isinstance(denom, PolyElement)
        assert numer.ring == ring.numer_ring
        assert denom.ring == ring.numer_ring
        assert denom.is_term

        [(denom_monom, denom_coeff)] = denom.iterterms()
        assert ring.domain.is_one(denom_coeff)

        min_degrees: tuple[int, ...] = tuple(map(min, zip(*numer.itermonoms())))
        for d, m in zip(denom_monom, min_degrees):
            assert d >= 0 and m >= 0, "Negative exponents in numer or denom"
            if check_cancelled:
                assert d == 0 or m == 0, "Uncancelled exponents between numer and denom"

    def new(self,
            numer: PolyElement,
            denom: Optional[PolyElement] = None
            ) -> 'LaurentPolyElement':
        """Construct a new Laurent polynomial in the same ring as ``self``."""
        ring = self.ring
        if denom is None:
            denom = ring.numer_ring.one
        numer, denom = self._normalize_poly(ring, numer, denom)
        return self.__class__(ring, numer, denom)

    @classmethod
    def _normalize_poly(cls,
                        ring: LaurentPolyRing,
                        numer: PolyElement,
                        denom: PolyElement
                        ) -> tuple[PolyElement, PolyElement]:
        """Normalize a poly/denom representation of a Laurent polynomial.

        The monomials are required to have non-negative exponents.
        """
        cls._check(ring, numer, denom, check_cancelled=False)

        [denom_monom] = denom.itermonoms()

        numer_ring = ring.numer_ring
        if denom_monom == numer_ring.zero_monom:
            return numer, denom

        monomial_gcd = numer_ring.monomial_gcd
        monomial_ldiv = numer_ring.monomial_ldiv
        domain = numer_ring.domain

        monom_gcd = denom_monom

        for monom in numer.itermonoms():
            monom_gcd = monomial_gcd(monom_gcd, monom)
            if monom_gcd == numer_ring.zero_monom:
                return numer, denom

        denom_monom_new = monomial_ldiv(denom_monom, monom_gcd)
        numer_new = {monomial_ldiv(m, monom_gcd): c for m, c in numer.items()}

        denom_poly = numer_ring.term_new(denom_monom_new, domain.one)
        numer_poly = numer_ring.from_dict(numer_new)

        return numer_poly, denom_poly

    def convert_to_ring(self, element: Any) -> Optional['LaurentPolyElement']:
        """Convert ``element`` to the ring of ``self``.

        Returns ``None`` if the conversion is not possible.
        """
        if isinstance(element, LaurentPolyElement) and self.ring == element.ring:
            return element
        elif isinstance(element, PolyElement) and self.ring.numer_ring == element.ring:
            return self.new(element)
        else:
            try:
                return self.ring.ground_new(element)
            except CoercionFailed:
                return None

    def _terms_num_den(self) -> list[tuple[PolyElement, PolyElement]]:
        denom = self.denom
        term_new = self.ring.numer_ring.term_new
        terms = [term_new(m, c) for m, c in self.numer.terms()]
        terms_num_den = []
        for term in terms:
            _, num, den = term._gcd_monom(denom)
            terms_num_den.append((num, den))
        return terms_num_den

    def as_expr(self, fraction=True) -> Expr:
        """Convert a Laurent polynomial to a SymPy expression."""
        if fraction:
            return self.as_expr_fraction()
        else:
            return self.as_expr_add()

    def as_expr_add(self) -> Expr:
        """Convert a Laurent polynomial to a SymPy expression."""
        terms = []
        for num, den in self._terms_num_den():
            if den == self.ring.numer_ring.one:
                terms.append(num.as_expr())
            else:
                terms.append(num.as_expr() / den.as_expr())
        return Add(*terms)

    def as_expr_fraction(self) -> Expr:
        """Convert a Laurent polynomial to a SymPy expression."""
        return self.numer.as_expr() / self.denom.as_expr()

    def __str__(self):
        terms_str = []
        for num, den in self._terms_num_den():
            if den == self.ring.numer_ring.one:
                terms_str.append(str(num))
            else:
                terms_str.append(f'{num}/{den}')

        if not terms_str:
            return '0'

        return ' + '.join(terms_str)

    @property
    def is_term(self):
        return self.numer.is_term

    @property
    def is_zero(self):
        return not self

    @property
    def is_one(self):
        return self.numer.is_one and self.denom.is_one

    @property
    def is_ground(self):
        return self.numer.is_ground and self.denom.is_ground

    @property
    def LC(self):
        return self.numer.LC

    def __bool__(self):
        return bool(self.numer)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.new(-self.numer, self.denom)

    def __add__(self, other: PolyElement) -> 'LaurentPolyElement':
        other_ring = self.convert_to_ring(other)
        if other_ring is None:
            return NotImplemented
        return self._add(other_ring)

    def __radd__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other._add(self)

    def _add(self, other):
        numer = self.numer * other.denom + self.denom * other.numer
        denom = self.denom * other.denom
        return self.new(numer, denom)

    def __sub__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return self._sub(other)

    def __rsub__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other._sub(self)

    def _sub(self, other):
        numer = self.numer * other.denom - self.denom * other.numer
        denom = self.denom * other.denom
        return self.new(numer, denom)

    def __mul__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return self._mul(other)

    def __rmul__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other._mul(self)

    def _mul(self, other):
        numer = self.numer * other.numer
        denom = self.denom * other.denom
        return self.new(numer, denom)

    def __pow__(self, n: int):
        if not isinstance(n, int):
            return NotImplemented

        if n == 0:
            return self.ring.one

        numer, denom = self.numer, self.denom

        if n < 0:
            if not self:
                raise ZeroDivisionError("Inversion of zero")
            elif not numer.is_term:
                raise NotImplementedError("Inversion of non-term")
            numer, denom = denom, numer
            n = -n

        return self.new(numer ** n, denom ** n)

    def __rtruediv__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other._div_term(self)

    def __truediv__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return self._div_term(other)

    def _div_term(self, other):
        if other.is_zero:
            raise ZeroDivisionError("Division by zero")
        if not other.is_term:
            raise NotImplementedError("Division by non-term")

        [(monom, coeff)] = other.numer.items()
        numer_ring = self.ring.numer_ring
        monom = numer_ring.term_new(monom, numer_ring.domain.one)

        # XXX: Use exquo_ground when available
        numer = (self.numer * other.denom).quo_ground(coeff)
        denom = self.denom * monom

        return self.new(numer, denom)

    def __mod__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return self.rem(other)

    def __rmod__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other.rem(self)

    def __floordiv__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return self.quo(other)

    def __rfloordiv__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other.quo(self)

    def __divmod__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return self.div(other)

    def __rdivmod__(self, other):
        other = self.convert_to_ring(other)
        if other is None:
            return NotImplemented
        return other.div(self)

    def div(self, other):
        """Quotient and remainder division of Laurent polynomials."""
        # self = other * quo + rem
        # p1/m1 = p2/m2 * q + r
        # p1 = p2*(q*m1/m2) + r*m1
        # q*m1/m2 = p1 // p2
        # r*m1 = p1 % p2
        # q = (p1 // p2) * (m2/m1)
        # r = (p1 % p2) / m1
        # deg(r*m1) < deg(p2)
        quo, rem = self.numer.div(other.numer)
        quo_laurent = self.new(quo * other.denom, self.denom)
        rem_laurent = self.new(rem, self.denom)
        return quo_laurent, rem_laurent

    def quo(self, other):
        return self.div(other)[0]

    def rem(self, other):
        return self.div(other)[1]

    def exquo(self, other):
        quo_numer = self.numer.exquo(other.numer)
        numer = quo_numer * other.denom
        denom = self.denom
        return self.new(numer, denom)

    def gcd(self, other):
        gcd_numer = self.numer.gcd(other.numer)
        gcd_denom = self.denom.gcd(other.denom)
        return self.new(gcd_numer, gcd_denom)

    def _factor_list(self):
        """Compute a factorization of ``self``."""
        # XXX: Should negative multiplicity should be used for the denominator
        # rather than positive powers of a reciprocal?
        #
        # We will leave this as a private method for now...
        coeff, numer_factors = self.numer.factor_list()
        _, denom_factors = self.denom.factor_list()

        one = self.ring.numer_ring.one
        factors = []
        for factor, exp in numer_factors:
            factor_laurent = self.new(factor, one)
            factors.append((factor_laurent, exp))
        for factor, exp in denom_factors:
            factor_laurent = self.new(one, factor)
            factors.append((factor_laurent, exp))

        return coeff, factors
