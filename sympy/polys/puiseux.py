"""
Puiseux rings. These are used by the ring_series module to represented
truncated Puiseux series. Elements of a Puiseux ring are like polynomials
except that the exponents can be negative or rational rather than just
non-negative integers.
"""

# Previously the ring_series module used PolyElement to represent Puiseux
# series. This is problematic because it means that PolyElement has to support
# negative and non-integer exponents which most polynomial representations do
# not support. This module provides an implementation of a ring for Puiseux
# series that can be used by ring_series without breaking the basic invariants
# of polynomial rings.
#
# Ideally there would be more of a proper series type that can keep track of
# not just the leading terms of a truncated series but also the precision
# of the series. For now the rings here are just introduced to keep the
# interface that ring_series was using before.

from __future__ import annotations

from typing import TYPE_CHECKING, Generic, Protocol, TypeVar

from math import gcd, lcm

from sympy.external.gmpy import MPQ
from sympy.polys.domains.domain import Er
from sympy.polys.domains import QQ
from sympy.polys.rings import PolyRing, PolyElement
from sympy.core.add import Add
from sympy.core.mul import Mul


if TYPE_CHECKING:
    from typing import Unpack
    from sympy.core.expr import Expr
    from sympy.polys.domains import Domain
    from collections.abc import Iterable, Iterator


K = TypeVar("K", covariant=True)
V = TypeVar("V", covariant=True)


class Map(Protocol[K, V]):
    """A dict mapping from keys to values."""

    def items(self) -> Iterable[tuple[K, V]]: ...
    def __iter__(self) -> Iterator[K]: ...


MonI = tuple[int, ...]
MonQ = tuple[MPQ, ...]


def puiseux_ring(
    symbols: str | list[Expr],
    domain: Domain[Er],
    /,
) -> tuple[PuiseuxRing[Er], Unpack[tuple[PuiseuxPoly[Er], ...]]]:
    """Construct a Puiseux ring.

    This function constructs a Puiseux ring with the given symbols and domain.

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.puiseux import puiseux_ring
    >>> R, x, y = puiseux_ring('x y', QQ)
    >>> R
    PuiseuxRing((x, y), QQ)
    >>> p = 5*x**QQ(1,2) + 7/y
    >>> p
    7*y**(-1) + 5*x**(1/2)
    """
    ring = PuiseuxRing(symbols, domain)
    return (ring,) + ring.gens  # type: ignore


class PuiseuxRing(Generic[Er]):
    """Ring of Puiseux polynomials.

    A Puiseux polynomial is a truncated Puiseux series. The exponents of the
    monomials can be negative or rational numbers. This ring is used by the
    ring_series module:

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.puiseux import puiseux_ring
    >>> from sympy.polys.ring_series import rs_exp, rs_nth_root
    >>> ring, x, y = puiseux_ring('x y', QQ)
    >>> f = x**2 + y**3
    >>> f
    y**3 + x**2
    >>> f.diff(x)
    2*x
    >>> rs_exp(x, x, 5)
    1 + x + 1/2*x**2 + 1/6*x**3 + 1/24*x**4

    Importantly the Puiseux ring can represent truncated series with negative
    and fractional exponents:

    >>> f = 1/x + 1/y**2
    >>> f
    x**(-1) + y**(-2)
    >>> f.diff(x)
    -1*x**(-2)

    >>> rs_nth_root(8*x + x**2 + x**3, 3, x, 5)
    2*x**(1/3) + 1/12*x**(4/3) + 23/288*x**(7/3) + -139/20736*x**(10/3)

    See Also
    ========

    sympy.polys.ring_series.rs_series
    PuiseuxPoly
    """

    def __init__(self, symbols: str | list[Expr], domain: Domain[Er]):
        poly_ring = PolyRing(symbols, domain)

        domain = poly_ring.domain
        ngens = poly_ring.ngens

        self.poly_ring = poly_ring
        self.domain = domain

        self.symbols = poly_ring.symbols
        self.gens = tuple([self.from_poly(g) for g in poly_ring.gens])
        self.ngens = ngens

        self.zero = self.from_poly(poly_ring.zero)
        self.one = self.from_poly(poly_ring.one)

        self.zero_monom = poly_ring.zero_monom
        self.monomial_mul = poly_ring.monomial_mul

    def __repr__(self) -> str:
        return f"PuiseuxRing({self.symbols}, {self.domain})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PuiseuxRing):
            return NotImplemented
        return self.symbols == other.symbols and self.domain == other.domain

    def from_poly(self, poly: PolyElement[Er]) -> PuiseuxPoly[Er]:
        """Create a Puiseux polynomial from a polynomial.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R1, x1 = ring('x', QQ)
        >>> R2, x2 = puiseux_ring('x', QQ)
        >>> R2.from_poly(x1**2)
        x**2
        """
        return PuiseuxPoly(poly, self)

    def from_dict(
        self, terms: Map[tuple[MPQ | int, ...], Er | Expr | int]
    ) -> PuiseuxPoly[Er]:
        """Create a Puiseux polynomial from a dictionary of terms.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x = puiseux_ring('x', QQ)
        >>> R.from_dict({(QQ(1,2),): QQ(3)})
        3*x**(1/2)
        """
        terms_q = {_to_monq(m): coeff for m, coeff in terms.items()}
        return PuiseuxPoly.from_dict(terms_q, self)

    def from_int(self, n: int) -> PuiseuxPoly[Er]:
        """Create a Puiseux polynomial from an integer.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x = puiseux_ring('x', QQ)
        >>> R.from_int(3)
        3
        """
        return self.from_poly(self.poly_ring(n))

    def domain_new(self, arg: Expr | Er | int) -> Er:
        """Create a new element of the domain.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x = puiseux_ring('x', QQ)
        >>> R.domain_new(3)
        3
        >>> QQ.of_type(_)
        True
        """
        return self.poly_ring.domain_new(arg)

    def ground_new(self, arg: Expr | Er | int) -> PuiseuxPoly[Er]:
        """Create a new element from a ground element.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.puiseux import puiseux_ring, PuiseuxPoly
        >>> R, x = puiseux_ring('x', QQ)
        >>> R.ground_new(3)
        3
        >>> isinstance(_, PuiseuxPoly)
        True
        """
        return self.from_poly(self.poly_ring.ground_new(arg))

    def __call__(
        self,
        arg: Expr
        | PolyElement[Er]
        | Er
        | int
        | Map[tuple[MPQ | int, ...], Er | Expr | int],
    ) -> PuiseuxPoly[Er]:
        """Coerce an element into the ring.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x = puiseux_ring('x', QQ)
        >>> R(3)
        3
        >>> R({(QQ(1,2),): QQ(3)})
        3*x**(1/2)
        """
        if isinstance(arg, dict):
            return self.from_dict(arg)
        else:
            return self.from_poly(self.poly_ring(arg))

    def index(self, x: PuiseuxPoly[Er]) -> int:
        """Return the index of a generator.

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x, y = puiseux_ring('x y', QQ)
        >>> R.index(x)
        0
        >>> R.index(y)
        1
        """
        return self.gens.index(x)


def _to_monq(monom: tuple[MPQ | int, ...]) -> MonQ:
    """Convert a tuple of MPQ or int to a tuple of MPQ."""
    return tuple(MPQ(mi) for mi in monom)


def _div_poly_monom(poly: PolyElement[Er], monom: MonI) -> PolyElement[Er]:
    ring = poly.ring
    div = ring.monomial_div
    return ring.from_dict({div(m, monom): c for m, c in poly.terms()})


def _mul_poly_monom(poly: PolyElement[Er], monom: MonI) -> PolyElement[Er]:
    ring = poly.ring
    mul = ring.monomial_mul
    return ring.from_dict({mul(m, monom): c for m, c in poly.terms()})


def _div_monom(monom: Iterable[int], div: Iterable[int]) -> MonI:
    return tuple(mi - di for mi, di in zip(monom, div))


class PuiseuxPoly(Generic[Er]):
    """Puiseux polynomial. Represents a truncated Puiseux series.

    See the :class:`PuiseuxRing` class for more information.

    >>> from sympy import QQ
    >>> from sympy.polys.puiseux import puiseux_ring
    >>> R, x, y = puiseux_ring('x, y', QQ)
    >>> p = 5*x**2 + 7*y**3
    >>> p
    7*y**3 + 5*x**2

    The internal representation of a Puiseux polynomial wraps a normal
    polynomial. To support negative powers the polynomial is considered to be
    divided by a monomial.

    >>> p2 = 1/x + 1/y**2
    >>> p2.monom # x*y**2
    (1, 2)
    >>> p2.poly
    x + y**2
    >>> (y**2 + x) / (x*y**2) == p2
    True

    To support fractional powers the polynomial is considered to be a function
    of ``x**(1/nx), y**(1/ny), ...``. The representation keeps track of a
    monomial and a list of exponent denominators so that the polynomial can be
    used to represent both negative and fractional powers.

    >>> p3 = x**QQ(1,2) + y**QQ(2,3)
    >>> p3.ns
    (2, 3)
    >>> p3.poly
    x + y**2

    See Also
    ========

    sympy.polys.puiseux.PuiseuxRing
    sympy.polys.rings.PolyElement
    """

    ring: PuiseuxRing[Er]
    poly: PolyElement[Er]
    monom: MonI | None
    ns: MonI | None

    def __new__(cls, poly: PolyElement[Er], ring: PuiseuxRing[Er]) -> PuiseuxPoly[Er]:
        return cls._new(ring, poly, None, None)

    @classmethod
    def _new(
        cls,
        ring: PuiseuxRing[Er],
        poly: PolyElement[Er],
        monom: MonI | None,
        ns: MonI | None,
    ) -> PuiseuxPoly[Er]:
        poly, monom, ns = cls._normalize(poly, monom, ns)
        return cls._new_raw(ring, poly, monom, ns)

    @classmethod
    def _new_raw(
        cls,
        ring: PuiseuxRing[Er],
        poly: PolyElement[Er],
        monom: MonI | None,
        ns: MonI | None,
    ) -> PuiseuxPoly[Er]:
        obj = object.__new__(cls)
        obj.ring = ring
        obj.poly = poly
        obj.monom = monom
        obj.ns = ns
        return obj

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PuiseuxPoly):
            return (
                self.poly == other.poly
                and self.monom == other.monom
                and self.ns == other.ns
            )
        elif self.monom is None and self.ns is None:
            return self.poly.__eq__(other)
        else:
            return NotImplemented

    @classmethod
    def _normalize(
        cls,
        poly: PolyElement[Er],
        monom: MonI | None,
        ns: MonI | None,
    ) -> tuple[PolyElement[Er], MonI | None, MonI | None]:
        if monom is None and ns is None:
            return poly, None, None

        if monom is not None:
            degs: MonI = tuple([max(d, 0) for d in poly.tail_degrees()])  # type: ignore
            if all(di >= mi for di, mi in zip(degs, monom)):
                poly = _div_poly_monom(poly, monom)
                monom = None
            elif any(degs):
                poly = _div_poly_monom(poly, degs)
                monom = _div_monom(monom, degs)

        if ns is not None:
            factors_d, [poly_d] = poly.deflate()
            degrees = poly.degrees()
            monom_d = monom if monom is not None else [0] * len(degrees)
            ns_new = []
            monom_new = []
            inflations = []
            for fi, ni, di, mi in zip(factors_d, ns, degrees, monom_d):
                if di == 0:
                    g = gcd(ni, mi)
                else:
                    g = gcd(fi, ni, mi)
                ns_new.append(ni // g)
                monom_new.append(mi // g)
                inflations.append(fi // g)

            if any(infl > 1 for infl in inflations):
                poly_d = poly_d.inflate(inflations)

            poly = poly_d

            if monom is not None:
                monom = tuple(monom_new)

            if all(n == 1 for n in ns_new):
                ns = None
            else:
                ns = tuple(ns_new)

        return poly, monom, ns

    @classmethod
    def _monom_fromint(
        cls,
        monom: MonI,
        dmonom: MonI | None,
        ns: MonI | None,
    ) -> MonQ:
        if dmonom is not None and ns is not None:
            return tuple(QQ(mi - di, ni) for mi, di, ni in zip(monom, dmonom, ns))
        elif dmonom is not None:
            return tuple(QQ(mi - di) for mi, di in zip(monom, dmonom))
        elif ns is not None:
            return tuple(QQ(mi, ni) for mi, ni in zip(monom, ns))
        else:
            return tuple(QQ(mi) for mi in monom)

    @classmethod
    def _monom_toint(
        cls,
        monom: MonQ,
        dmonom: MonI | None,
        ns: MonI | None,
    ) -> MonI:
        if dmonom is not None and ns is not None:
            return tuple(
                int((mi * ni).numerator + di) for mi, di, ni in zip(monom, dmonom, ns)
            )
        elif dmonom is not None:
            return tuple(int(mi.numerator + di) for mi, di in zip(monom, dmonom))
        elif ns is not None:
            return tuple(int((mi * ni).numerator) for mi, ni in zip(monom, ns))
        else:
            return tuple(int(mi.numerator) for mi in monom)

    def itermonoms(self) -> Iterator[MonQ]:
        """Iterate over the monomials of a Puiseux polynomial.

        >>> from sympy import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x, y = puiseux_ring('x, y', QQ)
        >>> p = 5*x**2 + 7*y**3
        >>> list(p.itermonoms())
        [(2, 0), (0, 3)]
        >>> p[(2, 0)]
        5
        """
        monom, ns = self.monom, self.ns
        for m in self.poly.itermonoms():
            yield self._monom_fromint(m, monom, ns)

    def monoms(self) -> list[MonQ]:
        """Return a list of the monomials of a Puiseux polynomial."""
        return list(self.itermonoms())

    def __iter__(self) -> Iterator[MonQ]:
        return self.itermonoms()

    def __getitem__(self, monom: MonQ) -> Er:
        monomq = self._monom_toint(monom, self.monom, self.ns)
        return self.poly[monomq]

    def __len__(self) -> int:
        return len(self.poly)

    def iterterms(self) -> Iterator[tuple[MonQ, Er]]:
        """Iterate over the terms of a Puiseux polynomial.

        >>> from sympy import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x, y = puiseux_ring('x, y', QQ)
        >>> p = 5*x**2 + 7*y**3
        >>> list(p.iterterms())
        [((2, 0), 5), ((0, 3), 7)]
        """
        monom, ns = self.monom, self.ns
        for m, coeff in self.poly.iterterms():
            mq = self._monom_fromint(m, monom, ns)
            yield mq, coeff

    def terms(self) -> list[tuple[MonQ, Er]]:
        """Return a list of the terms of a Puiseux polynomial."""
        return list(self.iterterms())

    @property
    def is_term(self) -> bool:
        """Return True if the Puiseux polynomial is a single term."""
        return self.poly.is_term

    def to_dict(self) -> dict[MonQ, Er]:
        """Return a dictionary representation of a Puiseux polynomial."""
        return dict(self.iterterms())

    @classmethod
    def from_dict(
        cls,
        terms: Map[MonQ, Er | Expr | int],
        ring: PuiseuxRing[Er],
    ) -> PuiseuxPoly[Er]:
        """Create a Puiseux polynomial from a dictionary of terms.

        >>> from sympy import QQ
        >>> from sympy.polys.puiseux import puiseux_ring, PuiseuxPoly
        >>> R, x = puiseux_ring('x', QQ)
        >>> PuiseuxPoly.from_dict({(QQ(1,2),): QQ(3)}, R)
        3*x**(1/2)
        >>> R.from_dict({(QQ(1,2),): QQ(3)})
        3*x**(1/2)
        """
        ns = [1] * ring.ngens
        mon: list[MPQ | int] = [0] * ring.ngens
        for mo in terms:
            ns = [lcm(n, m.denominator) for n, m in zip(ns, mo)]
            mon = [min(m, n) for m, n in zip(mo, mon)]

        if not any(mon):
            monom = None
        else:
            monom = tuple(-int((m * n).numerator) for m, n in zip(mon, ns))

        if all(n == 1 for n in ns):
            ns_final = None
        else:
            ns_final = tuple(ns)

        def conv(m: MonQ) -> MonI:
            return cls._monom_toint(m, monom, ns_final)

        terms_p = {conv(m): coeff for m, coeff in terms.items()}

        poly = ring.poly_ring.from_dict(terms_p)

        return cls._new(ring, poly, monom, ns_final)

    def as_expr(self) -> Expr:
        """Convert a Puiseux polynomial to :class:`~sympy.core.expr.Expr`.

        >>> from sympy import QQ, Expr
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x = puiseux_ring('x', QQ)
        >>> p = 5*x**2 + 7*x**3
        >>> p.as_expr()
        7*x**3 + 5*x**2
        >>> isinstance(_, Expr)
        True
        """
        ring = self.ring
        dom = ring.domain
        symbols = ring.symbols
        terms = []
        for monom, coeff in self.iterterms():
            coeff_expr = dom.to_sympy(coeff)
            monoms_expr = []
            for i, m in enumerate(monom):
                monoms_expr.append(symbols[i] ** m)
            terms.append(Mul(coeff_expr, *monoms_expr))
        return Add(*terms)

    def __repr__(self) -> str:
        def format_power(base: str, exp: MPQ) -> str:
            if exp == 1:
                return base
            elif exp >= 0 and int(exp) == exp:
                return f"{base}**{exp}"
            else:
                return f"{base}**({exp})"

        ring = self.ring
        dom = ring.domain

        syms = [str(s) for s in ring.symbols]
        terms_str = []
        for monom, coeff in sorted(self.terms()):
            monom_str = "*".join(format_power(s, e) for s, e in zip(syms, monom) if e)
            if coeff == dom.one:
                if monom_str:
                    terms_str.append(monom_str)
                else:
                    terms_str.append("1")
            elif not monom_str:
                terms_str.append(str(coeff))
            else:
                terms_str.append(f"{coeff}*{monom_str}")

        return " + ".join(terms_str)

    def _unify(
        self, other: PuiseuxPoly[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], MonI | None, MonI | None]:
        """Bring two Puiseux polynomials to a common monom and ns."""
        poly1, monom1, ns1 = self.poly, self.monom, self.ns
        poly2, monom2, ns2 = other.poly, other.monom, other.ns

        if monom1 == monom2 and ns1 == ns2:
            return poly1, poly2, monom1, ns1

        if ns1 == ns2:
            ns = ns1
        elif ns1 is not None and ns2 is not None:
            ns = tuple(lcm(n1, n2) for n1, n2 in zip(ns1, ns2))
            f1 = [n // n1 for n, n1 in zip(ns, ns1)]
            f2 = [n // n2 for n, n2 in zip(ns, ns2)]
            poly1 = poly1.inflate(f1)
            poly2 = poly2.inflate(f2)
            if monom1 is not None:
                monom1 = tuple(m * f for m, f in zip(monom1, f1))
            if monom2 is not None:
                monom2 = tuple(m * f for m, f in zip(monom2, f2))
        elif ns2 is not None:
            ns = ns2
            poly1 = poly1.inflate(ns)
            if monom1 is not None:
                monom1 = tuple(m * n for m, n in zip(monom1, ns))
        elif ns1 is not None:
            ns = ns1
            poly2 = poly2.inflate(ns)
            if monom2 is not None:
                monom2 = tuple(m * n for m, n in zip(monom2, ns))
        else:
            assert False

        if monom1 == monom2:
            monom = monom1
        elif monom1 is not None and monom2 is not None:
            monom = tuple(max(m1, m2) for m1, m2 in zip(monom1, monom2))
            poly1 = _mul_poly_monom(poly1, _div_monom(monom, monom1))
            poly2 = _mul_poly_monom(poly2, _div_monom(monom, monom2))
        elif monom2 is not None:
            monom = monom2
            poly1 = _mul_poly_monom(poly1, monom2)
        elif monom1 is not None:
            monom = monom1
            poly2 = _mul_poly_monom(poly2, monom1)
        else:
            assert False

        return poly1, poly2, monom, ns

    def __pos__(self) -> PuiseuxPoly[Er]:
        return self

    def __neg__(self) -> PuiseuxPoly[Er]:
        return self._new_raw(self.ring, -self.poly, self.monom, self.ns)

    def __add__(self, other: PuiseuxPoly[Er] | Er | int) -> PuiseuxPoly[Er]:
        if isinstance(other, PuiseuxPoly):
            if self.ring != other.ring:
                raise ValueError("Cannot add Puiseux polynomials from different rings")
            return self._add(other)
        domain = self.ring.domain
        if isinstance(other, int):
            return self._add_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            return NotImplemented

    def __radd__(self, other: Er | int) -> PuiseuxPoly[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._add_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            return NotImplemented

    def __sub__(self, other: PuiseuxPoly[Er] | Er | int) -> PuiseuxPoly[Er]:
        if isinstance(other, PuiseuxPoly):
            if self.ring != other.ring:
                raise ValueError(
                    "Cannot subtract Puiseux polynomials from different rings"
                )
            return self._sub(other)
        domain = self.ring.domain
        if isinstance(other, int):
            return self._sub_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._sub_ground(other)
        else:
            return NotImplemented

    def __rsub__(self, other: Er | int) -> PuiseuxPoly[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._rsub_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._rsub_ground(other)
        else:
            return NotImplemented

    def __mul__(self, other: PuiseuxPoly[Er] | Er | int) -> PuiseuxPoly[Er]:
        if isinstance(other, PuiseuxPoly):
            if self.ring != other.ring:
                raise ValueError(
                    "Cannot multiply Puiseux polynomials from different rings"
                )
            return self._mul(other)
        domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            return NotImplemented

    def __rmul__(self, other: Er | int) -> PuiseuxPoly[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            return NotImplemented

    def __pow__(self, other: int | MPQ) -> PuiseuxPoly[Er]:
        if isinstance(other, int):
            if other >= 0:
                return self._pow_pint(other)
            else:
                return self._pow_nint(-other)
        elif QQ.of_type(other):
            return self._pow_rational(other)
        else:
            return NotImplemented

    def __truediv__(self, other: PuiseuxPoly[Er] | Er | int) -> PuiseuxPoly[Er]:
        if isinstance(other, PuiseuxPoly):
            if self.ring != other.ring:
                raise ValueError(
                    "Cannot divide Puiseux polynomials from different rings"
                )
            return self._mul(other._inv())
        domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(domain.convert_from(QQ(1, other), QQ))
        elif domain.of_type(other):
            return self._div_ground(other)
        else:
            return NotImplemented

    def __rtruediv__(self, other: Er | int) -> PuiseuxPoly[Er]:
        if isinstance(other, int):
            return self._inv()._mul_ground(self.ring.domain.convert_from(QQ(other), QQ))
        elif self.ring.domain.of_type(other):
            return self._inv()._mul_ground(other)
        else:
            return NotImplemented

    def _add(self, other: PuiseuxPoly[Er]) -> PuiseuxPoly[Er]:
        poly1, poly2, monom, ns = self._unify(other)
        return self._new(self.ring, poly1 + poly2, monom, ns)

    def _add_ground(self, ground: Er) -> PuiseuxPoly[Er]:
        return self._add(self.ring.ground_new(ground))

    def _sub(self, other: PuiseuxPoly[Er]) -> PuiseuxPoly[Er]:
        poly1, poly2, monom, ns = self._unify(other)
        return self._new(self.ring, poly1 - poly2, monom, ns)

    def _sub_ground(self, ground: Er) -> PuiseuxPoly[Er]:
        return self._sub(self.ring.ground_new(ground))

    def _rsub_ground(self, ground: Er) -> PuiseuxPoly[Er]:
        return self.ring.ground_new(ground)._sub(self)

    def _mul(self, other: PuiseuxPoly[Er]) -> PuiseuxPoly[Er]:
        poly1, poly2, monom, ns = self._unify(other)
        if monom is not None:
            monom = tuple(2 * e for e in monom)
        return self._new(self.ring, poly1 * poly2, monom, ns)

    def _mul_ground(self, ground: Er) -> PuiseuxPoly[Er]:
        return self._new_raw(self.ring, self.poly * ground, self.monom, self.ns)

    def _div_ground(self, ground: Er) -> PuiseuxPoly[Er]:
        return self._new_raw(self.ring, self.poly / ground, self.monom, self.ns)

    def _pow_pint(self, n: int) -> PuiseuxPoly[Er]:
        assert n >= 0
        monom = self.monom
        if monom is not None:
            monom = tuple(m * n for m in monom)
        return self._new(self.ring, self.poly**n, monom, self.ns)

    def _pow_nint(self, n: int) -> PuiseuxPoly[Er]:
        return self._inv()._pow_pint(n)

    def _pow_rational(self, n: MPQ) -> PuiseuxPoly[Er]:
        if not self.is_term:
            raise ValueError("Only monomials can be raised to a rational power")
        [(monom, coeff)] = self.terms()
        domain = self.ring.domain
        if not domain.is_one(coeff):
            raise ValueError("Only monomials can be raised to a rational power")
        monom = tuple(m * n for m in monom)
        return self.ring.from_dict({monom: domain.one})

    def _inv(self) -> PuiseuxPoly[Er]:
        if not self.is_term:
            raise ValueError("Only terms can be inverted")
        [(monom, coeff)] = self.terms()
        domain = self.ring.domain
        if not domain.is_Field and not domain.is_one(coeff):
            raise ValueError("Cannot invert non-unit coefficient")
        monom = tuple(-m for m in monom)
        coeff = 1 / coeff  # type: ignore
        return self.ring.from_dict({monom: coeff})

    def diff(self, x: PuiseuxPoly[Er]) -> PuiseuxPoly[Er]:
        """Differentiate a Puiseux polynomial with respect to a variable.

        >>> from sympy import QQ
        >>> from sympy.polys.puiseux import puiseux_ring
        >>> R, x, y = puiseux_ring('x, y', QQ)
        >>> p = 5*x**2 + 7*y**3
        >>> p.diff(x)
        10*x
        >>> p.diff(y)
        21*y**2
        """
        ring = self.ring
        i = ring.index(x)
        g: dict[MonQ, Er] = {}
        for expv, coeff in self.iterterms():
            n = expv[i]
            if n:
                e = list(expv)
                e[i] -= 1
                g[tuple(e)] = coeff * n # type: ignore
        return ring.from_dict(g)
