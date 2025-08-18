from sympy.polys.series.ring import PowerSeriesRing, PowerSeriesElement, TSeries
from sympy.polys.domains.ring import Ring
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domain import Er, Domain
from sympy.core.expr import Expr
from sympy.utilities import public
from typing import Type, TypeIs


@public
class Series(Ring[PowerSeriesElement[Er]], CompositeDomain):
    """A class for representing multivariate power series rings."""

    is_PowerSeriesRing = is_Series = True

    has_assoc_Ring = True
    has_assoc_Field = False

    ring: PowerSeriesRing[Er]
    dtype: Type[PowerSeriesElement[Er]]
    gen: PowerSeriesElement[Er]
    symbol: Expr
    domain: Domain[Er]

    def __init__(self, domain: Domain[Er], symbol: str | Expr = "x", prec: int = 6):
        ring = PowerSeriesRing(domain, symbol, prec)

        self.ring = ring
        self.dtype = ring.dtype

        self.domain = ring.domain
        self.gen = ring.gen
        self.symbol = ring.symbol

    def new(self, element: TSeries | Expr | Er | int) -> PowerSeriesElement[Er]:  # type: ignore
        return self.ring.ring_new(element)

    def of_type(self, element) -> TypeIs[PowerSeriesElement[Er]]:
        """Check if ``a`` is of type ``dtype``."""
        return self.ring.is_element(element)

    @property
    def zero(self) -> PowerSeriesElement[Er]:  # type: ignore
        return self.ring.zero

    @property
    def one(self) -> PowerSeriesElement[Er]:  # type: ignore
        return self.ring.one

    @property
    def prec(self) -> int:
        return self.ring.prec

    def __repr__(self) -> str:
        return f"{self.domain}[[{self.symbol}], {self.prec}]"

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.ring, self.domain, self.symbol))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Series):
            return NotImplemented
        return self.ring == other.ring

    def is_unit(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns ``True`` if ``LC`` of series is a unit of ``self``"""
        if not a.is_ground:
            return False
        K = self.domain
        return K.is_unit(K.convert_from(a, self))

    def canonical_unit(self, a: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        u = self.domain.canonical_unit(a.leading_coefficient())
        return self.ring.ground_new(u)

    def to_sympy(self, a: PowerSeriesElement[Er]) -> Expr:
        """Convert `a` to a SymPy object."""
        return a.as_expr()

    def from_sympy(self, a: Expr) -> PowerSeriesElement[Er]:
        """Convert SymPy's expression to `dtype`."""
        return self.ring.from_expr(a)

    def from_ZZ(K1, a, K0):
        """Convert a Python `int` object to `dtype`."""
        return K1(K1.domain.convert(a, K0))

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`."""
        return K1(K1.domain.convert(a, K0))

    def from_QQ(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`."""
        return K1(K1.domain.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`."""
        return K1(K1.domain.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`."""
        return K1(K1.domain.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`."""
        return K1(K1.domain.convert(a, K0))

    def is_positive(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `LC(a)` is positive."""
        lc = a.leading_coefficient()
        return self.domain.is_positive(lc)

    def is_negative(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `LC(a)` is negative."""
        lc = a.leading_coefficient()
        return self.domain.is_negative(lc)

    def is_nonpositive(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `LC(a)` is non-positive."""
        lc = a.leading_coefficient()
        return self.domain.is_nonpositive(lc)

    def is_nonnegative(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `LC(a)` is non-negative."""
        lc = a.leading_coefficient()
        return self.domain.is_nonnegative(lc)
