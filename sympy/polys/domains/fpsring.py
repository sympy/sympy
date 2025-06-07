from __future__ import annotations

from collections.abc import Iterator
from typing import Any

from sympy import QQ
from sympy.core.expr import Expr
from sympy.polys.densearith import (dup_add, dup_add_ground, dup_mul,
                                    dup_mul_ground,dup_pow, dup_sub,
                                    dup_sub_ground)
from sympy.polys.domains import Domain
from sympy.polys.orderings import LexOrder, lex
from sympy.polys.polyclasses import DMP_Python
from sympy.polys.polyerrors import CoercionFailed, GeneratorsError
from sympy.polys.rings import _parse_symbols
from sympy.utilities import public


@public
class PowerSeriesRing:
    """A Domain for Univariate Power Series Ring."""

    is_PowerSeriesRing: bool = True
    is_PowerSeries: bool = True
    has_assoc_Ring: bool = True
    has_assoc_Field: bool = False
    DEFAULT_PRECISION: int = 6

    domain: Domain
    symbol: str | tuple[Expr, ...] | list[Expr]
    prec: int
    order: LexOrder

    def __init__(
        self,
        domain: Domain,
        symbol: str | list[Expr] | tuple[Expr, ...],
        prec: int | None = None,
        ) -> None:
        """Initialize a PowerSeriesRing.

        Parameters
        ==========
        domain : Domain
            The domain of the coefficients of the power series
        symbol : str or list of Expr or tuple of Expr
            The variables used in the power series
        prec : int, optional
            The precision (maximum degree) used for power series operations.
            If None, uses DEFAULT_PRECISION.
        """

        symbol = tuple(_parse_symbols(symbol))
        if len(symbol) != 1:
            raise GeneratorsError("Univariate Power series ring must have exactly one symbol.")

        if prec is None:
            prec = self.DEFAULT_PRECISION
        elif prec < 1:
            raise ValueError("Precision must be a positive integer.")


        self.domain = domain
        self.symbol = symbol
        self.prec = prec
        self.order = lex


    @property
    def zero(self) -> PowerSeriesElement:
        """Return the zero element of the power series ring."""
        return self.from_list([0])

    @property
    def one(self) -> PowerSeriesElement:
        """Return the one element of the power series ring."""
        return self.from_list([1])

    @property
    def ring(self) -> str:
        """Representation of working ring."""
        return f"Power Series Ring in {self.symbol[0]} over {self.domain} with {self.order} order"

    @property
    def gen(self) -> PowerSeriesElement:
        """Return the generator of the power series ring."""
        return self.from_list([1,0])


    def __repr__(self) -> str:
        return f'{self.domain}[[{", ".join(str(s) for s in self.symbol)}], {self.prec}]'

    def __str__(self) -> str:
        return f'{self.domain}[[{", ".join(str(s) for s in self.symbol)}], {self.prec}]'

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.domain, self.symbol, self.prec))

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, PowerSeriesRing)
            and self.domain == other.domain
            and self.symbol == other.symbol
            and self.prec == other.prec
        )

    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def ground_new(self, coeff: Any) -> PowerSeriesElement:
        """Create a power series element from a ground coefficient."""
        return self.term_new(coeff, 0)

    def term_new(self, coeff: Any, power: int) -> PowerSeriesElement:
        """Create a power series term less than precission of a ring."""
        if power < 0:
            raise ValueError("Power must be a non-negative integer.")
        elif power >= self.prec:
            raise ValueError(f"Power must be less than the precision {self.prec}.")
        elif not isinstance(power, int):
            raise ValueError("Power must be an integer.")

        new_series = list([coeff] + ([0] * (power-1)))
        return PowerSeriesElement(self, new_series)

    def from_list(self, series: list[Any]) -> PowerSeriesElement:
        """Create a power series element from a list of coefficients."""
        if not isinstance(series, list):
            raise TypeError("Series must be a list of coefficients.")
        return PowerSeriesElement(self, series)

    def order_term(self) -> PowerSeriesElement:
        """Return the order term of the power series."""
        return self.zero


class PowerSeriesElement:
    """Element of Univariate Power Series Ring. """

    ring: PowerSeriesRing
    domain: Domain
    poly: DMP_Python
    symbol: str | list[Expr] | tuple[Expr, ...]
    prec: int

    def __new__(cls,
                ring: PowerSeriesRing,
                series: list[Any]
                ) -> 'PowerSeriesElement':
        """Create a power series element.

        Parameters
        ==========
        ring : PowerSeriesRing
            The ring this element belongs to, defining the domain and variables
        series : list
            List of coefficients for the power series, where the index
            corresponds to the power of the variable
        """
        prec: int = ring.prec
        # Truncate to precision
        coeffs = series[-prec:] if len(series) > prec else series

        obj = object.__new__(cls)
        obj.ring = ring
        obj.poly = DMP_Python(coeffs, ring.domain, 0)
        return obj

    def _new(self, series: list[Any]) -> PowerSeriesElement:
        """Create a new PowerSeriesElement with the same domain and symbol."""
        return PowerSeriesElement(self.ring, series)

    def __hash__(self):
        """Hash function for the power series element."""
        return hash((self.__class__.__name__, self.ring, self.poly))

    def iterterms(self) -> Iterator[Any]:
        """Return an iterator over the terms of the power series."""
        return iter(reversed(self.poly.to_list()))

    def __repr__(self) -> str:
        """Representation string for the power series."""
        ring: PowerSeriesRing = self.ring
        prec: int = ring.prec
        sym: str | Expr = ring.symbol[0]
        domain: Domain = ring.domain
        poly_list = self.poly.to_list()

        if not poly_list or all(coeff == 0 for coeff in poly_list):
            return f"O({sym}**{prec})"

        series: str = ""
        first_term: bool = True
        terms: list[tuple[int, Any]] = []

        for pow, coeff in enumerate(self.iterterms()):
            if pow >= prec:
                break
            if coeff != 0:
                terms.append((pow, coeff))

        for pow, coeff in terms:
            negative: bool = domain.is_negative(coeff)
            sign: str = " - " if negative else " + "

            if first_term:
                if negative:
                    series += "-"
                first_term = False
            else:
                series += sign

            if negative:
                coeff = -coeff

            term: str
            if pow == 0:
                term = f"{coeff}"
            elif pow == 1:
                if coeff == 1:
                    term = f"{sym}"
                else:
                    term = f"{coeff}*{sym}"
            else:
                if coeff == 1:
                    term = f"{sym}**{pow}"
                else:
                    term = f"{coeff}*{sym}**{pow}"

            if negative and term.startswith("-"):
                term = term[1:]

            series += term

        return series + f" + O({sym}**{prec})"

    def __eq__(self, other: object) -> bool:
        """Equality check for two power series."""
        return (
            isinstance(other, PowerSeriesElement)
            and self.ring == other.ring
            and self.poly == other.poly
        )

    def __ne__(self, other: object) -> bool:
        """Inequality check for two power series."""
        return not self.__eq__(other)

    def __pos__(self) -> PowerSeriesElement:
        """Return the power series itself (unary plus)."""
        return self

    def __neg__(self) -> PowerSeriesElement:
        """Negation of the power series."""
        new_series = [-coeff for coeff in self.poly.to_list()]
        return self._new(new_series)

    def __add__(self, other: PowerSeriesElement) -> PowerSeriesElement:
        """Add another power series or scalar to this power series.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 5)
        >>> x = R.gen
        >>> f = 1 + 2*x + 3*x**2
        >>> g = 1 + 5*x**2
        >>> f + g
        2 + 2*x + 8*x**2 + O(x**5)
        """
        if isinstance(other, PowerSeriesElement):
            if self.ring != other.ring:
                raise ValueError("Cannot add PowerSeriesElement from different domains")
            return self._add(other)

        domain: Domain = self.ring.domain
        if isinstance(other, int):
            return self._add_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            raise NotImplementedError

    def __radd__(self, other: Any) -> PowerSeriesElement:
        """Add this power series to a scalar (right addition).

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 5)
        >>> x = R.gen
        >>> f = 1 + 2*x + 3*x**2
        >>> 3 + f
        4 + 2*x + 3*x**2 + O(x**5)
        """
        domain: Domain = self.ring.domain
        if isinstance(other, int):
            return self._add_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            raise NotImplementedError

    def __sub__(self, other: Any) -> PowerSeriesElement:
        """Subtract another power series or scalar from this power series.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 5)
        >>> x = R.gen
        >>> f = 1 + 2*x + 3*x**2
        >>> g = 1 + 5*x**2
        >>> f - g
        2*x - 2*x**2 + O(x**5)
        """
        if isinstance(other, PowerSeriesElement):
            if self.ring != other.ring:
                raise ValueError("Cannot subtract PowerSeriesElement from different domains")
            return self._sub(other)

        ring: PowerSeriesRing = self.ring
        domain: Domain = ring.domain

        try:
            if isinstance(other, int):
                other_coef = domain.convert_from(QQ(other), QQ)
            elif domain.of_type(other):
                other_coef = other
            else:
                raise NotImplementedError

            return self._sub_ground(other_coef)
        except (CoercionFailed, NotImplementedError):
            return NotImplemented

    def __rsub__(self, other: Any) -> PowerSeriesElement:
        """Subtract this power series from a scalar (right subtraction).

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 5)
        >>> x = R.gen
        >>> f = 1 + 2*x + 3*x**2
        >>> 3 - f
        2 - 2*x - 3*x**2 + O(x**5)
        """
        ring: PowerSeriesRing = self.ring
        domain: Domain = ring.domain

        try:
            if isinstance(other, int):
                other_coef = domain.convert_from(QQ(other), QQ)
            elif domain.of_type(other):
                other_coef = other
            else:
                raise NotImplementedError

            result = self.__neg__()

            return result._add_ground(other_coef)
        except (CoercionFailed, NotImplementedError):
            return NotImplemented

    def __mul__(self, other: int | PowerSeriesElement) -> PowerSeriesElement:
        """Multiply this power series by another power series or scalar.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 10)
        >>> f = 1 + 2*x + 3*x**2
        >>> g =  5 + 6*x + x**4
        >>> f * g
        5 + 16*x + 27*x**2 + 18*x**3 + x**4 + 2*x**5 + 3*x**6 + O(x**10)
        """
        if isinstance(other, PowerSeriesElement):
            if self.ring != other.ring:
                raise ValueError("Cannot multiply PowerSeriesElement from different domains")
            return self._mul(other)

        domain: Domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            raise NotImplementedError

    def __rmul__(self, other: int | Any) -> PowerSeriesElement:
        """Multiply a scalar by this power series (right multiplication).

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 10)
        >>> x = R.gen
        >>> f =  1 + 2*x + 3*x**2
        >>> 5 * f
        5 + 10*x + 15*x**2 + O(x**10)
        """
        domain: Domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(domain.convert_from(QQ(other), QQ))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            raise NotImplementedError

    def __pow__(self, n: int) -> PowerSeriesElement:
        """Raise this power series to an `n` power.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> R = PowerSeriesRing(QQ, 'x', 10)
        >>> x = R.gen
        >>> f =  1 + x + x**4
        >>> f**2
        1 + 2*x + x**2 + x**8 + O(x**10)
        """
        if not isinstance(n, int):
            raise ValueError("Exponent must be an integer.")
        if n < 0:
            raise ValueError("Negative exponent is not supported for power series.")
        return self._pow(n)

    def _add(self, other: PowerSeriesElement) -> PowerSeriesElement:
        new_series = dup_add(self.poly.to_list(),
                             other.poly.to_list(), self.ring.domain)
        return self._new(new_series)

    def _add_ground(self, other: PowerSeriesElement) -> PowerSeriesElement:
        new_series = dup_add_ground(self.poly.to_list(),
                                    other, self.ring.domain)
        return self._new(new_series)

    def _sub(self, other: PowerSeriesElement) -> PowerSeriesElement:
        new_series = dup_sub(self.poly.to_list(),
                             other.poly.to_list(), self.ring.domain)
        return self._new(new_series)

    def _sub_ground(self, other: PowerSeriesElement) -> PowerSeriesElement:
        new_series = dup_sub_ground(self.poly.to_list(),
                                    other, self.ring.domain)
        return self._new(new_series)

    def _mul(self, other: PowerSeriesElement) -> PowerSeriesElement:
        new_series = dup_mul(self.poly.to_list(),
                           other.poly.to_list(), self.ring.domain)
        return self._new(new_series)

    def _mul_ground(self, other: PowerSeriesElement) -> PowerSeriesElement:
        new_series = dup_mul_ground(self.poly.to_list(),
                                  other, self.ring.domain)
        return self._new(new_series)

    def _pow(self, n: int) -> PowerSeriesElement:
        new_series = dup_pow(self.poly.to_list(),
                           n, self.ring.domain)
        return self._new(new_series)

    def type(self) -> type:
        """Return the type of the power series element."""
        return PowerSeriesElement

    @property
    def is_ground(self) -> bool:
        """Check if the power series is a ground element (constant term only)."""
        return len(self.poly.to_list()) == 1 and self.poly.to_list()[0] != 0

    def truncate(self, n: int) -> PowerSeriesElement:
        """Truncate the power series to the first n terms.

        Examples
        ========
        >>> from sympy import QQ, symbols
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> x = symbols('x')
        >>> R = PowerSeriesRing(QQ, x, 5)
        >>> x = R.gen
        >>> f = 1 + 2*x + 3*x**2
        >>> f.truncate(2)
        1 + 2*x + O(x**5)
        """

        coeff: list[Any] = self.poly.to_list()
        prec: int = self.ring.prec

        if n < 0:
            raise ValueError("n must be a non-negative integer.")
        if n > prec:
            raise ValueError(f"Truncating n must be less than or equal to the precision {prec}.")
        new_series = coeff[-n:] if len(coeff) > n else coeff
        return self._new(new_series)

    def get_constant_term(self) -> Any:
        """Returns the coefficient of x^0 (the constant term) in the power series.

        Examples
        ========
        >>> from sympy import QQ, symbols
        >>> from sympy.polys.domains.fpsring import PowerSeriesRing
        >>> x = symbols('x')
        >>> R = PowerSeriesRing(QQ, x, 5)
        >>> x = R.gen
        >>> f = 1 + 2*x + 3*x**2
        >>> f.get_constant_term()
        1
        """

        coeff = self.poly.to_list()
        return coeff[-1] if coeff else self.ring.zero
