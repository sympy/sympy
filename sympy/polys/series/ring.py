from __future__ import annotations

from functools import reduce
from operator import add, mul

from sympy.core.expr import Expr
from sympy.core.symbol import Symbol
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.core.sympify import CantSympify, sympify
from sympy.polys.domains.domain import Er, Ef, Domain, DomainElement
from sympy.polys.densebasic import dup
from sympy.polys.domains.field import Field
from sympy.polys.polyerrors import GeneratorsError
from sympy.polys.series.base import PowerSeriesRingProto, PowerSeriesRingFieldProto
from sympy.polys.series.tring import TSeriesElement, _power_series_ring
from sympy.series.order import Order


from typing import Generic, overload, TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TypeIs


@overload
def power_series_ring(
    symbol: str, K: Field[Ef], prec: int = 6
) -> tuple[PowerSeriesRingField[Ef], PowerSeriesElement[Ef]]: ...


@overload
def power_series_ring(
    symbol: str, K: Domain[Er], prec: int = 6
) -> tuple[PowerSeriesRingRing[Er], PowerSeriesElement[Er]]: ...


def power_series_ring(
    symbol: str, K: Domain, prec: int = 6
) -> tuple[PowerSeriesRingRing | PowerSeriesRingField, PowerSeriesElement]:
    """
    Create a power series ring over the given domain.

    Parameters
    ==========

    symbol : str
        The symbol to use for the power series variable.

    K : Domain
        The ground domain for the power series ring. Must be ZZ or QQ.

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Returns
    =======

    ring : PowerSeriesRing
        The power series ring.

    generator : PowerSeriesElement
        The generator of the power series ring.

    Examples
    ========

    >>> from sympy import QQ
    >>> from sympy.polys.series import power_series_ring
    >>> R, x = power_series_ring("x", QQ)
    >>> R.sin(x + x**2)
    x + x**2 - 1/6*x**3 - 1/2*x**4 - 59/120*x**5 + O(x**6)
    >>> R.log(1 + 7*x**2)
    7*x**2 - 49/2*x**4 + O(x**6)

    See Also
    ========

    sympy.polys.series.ring.PowerSeriesRingRing
    sympy.polys.series.ring.PowerSeriesRingField
    sympy.polys.series.ring.PowerSeriesElement

    """

    if K.is_ZZ:
        ring = PowerSeriesRingRing(K, symbol, prec)
        return ring, ring.gen
    elif K.is_QQ:
        ring = PowerSeriesRingField(K, symbol, prec)
        return ring, ring.gen
    else:
        raise ValueError(f"Unsupported ground domain: {K}")


class PowerSeriesRingRing(Generic[Er]):
    """A class for representing Univariate Power Series Rings over a Ring."""

    ring: PowerSeriesRingProto[TSeriesElement[Er], Er]
    dtype: type[PowerSeriesElement[Er]]
    domain: Domain[Er]
    symbol: Expr
    prec: int

    def __init__(self, domain: Domain[Er], symbol: str | Expr = "x", prec: int = 6):
        if not isinstance(symbol, (str, Expr)):
            raise GeneratorsError("Symbol must be a string or a SymPy expression.")

        if isinstance(symbol, str):
            symbol = Symbol(symbol)

        self.ring = _power_series_ring(domain, prec)

        self.dtype = PowerSeriesElement
        self.domain = domain
        self.symbol = symbol
        self.prec = prec

        self.one = self.from_ground(self.domain.one)
        self.zero = self.from_ground(self.domain.zero)
        self.gen = self.from_element(self.ring.gen)

    def __repr__(self) -> str:
        return (
            f"Power Series Ring in {self.symbol} over {self.domain} of size {self.prec}"
        )

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PowerSeriesRingRing):
            return self.ring == other.ring and self.symbol == other.symbol
        else:
            return NotImplemented

    def __hash__(self) -> int:
        return hash((self.ring, self.domain, self.symbol, self.prec))

    def is_element(self, element: object) -> TypeIs[PowerSeriesElement[Er]]:
        """Check if element belongs to this ring."""
        return isinstance(element, PowerSeriesElement) and element.ring == self.ring

    def order_term(self):
        """
        Return the order term of the power series.

        Examples
        ========
        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x")
        >>> R.order_term()
        O(x**6)
        >>> R10 = PowerSeriesRingRing(ZZ, "x", 10)
        >>> R10.order_term()
        O(x**10)

        """
        o = self.ring([], self.prec)
        return self.from_element(o)

    def truncate(self, s: PowerSeriesElement[Er], prec: int) -> PowerSeriesElement[Er]:
        """Truncate the power series to the given precision."""
        R = self.ring
        return self.from_element(R.truncate(s.series, prec))

    def from_expr(self, expr: Expr) -> PowerSeriesElement[Er]:
        """
        Convert an expression to a power series element.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x", 5)
        >>> from sympy.abc import x
        >>> s = R.from_expr(x + x**2 + x**4 + x**7); s
        x + x**2 + x**4 + O(x**5)
        >>> type(s)
        <class 'sympy.polys.series.ring.PowerSeriesElement'>

        """
        return self._rebuild_expr(expr)

    def _rebuild_expr(self, expr) -> PowerSeriesElement[Er]:
        """Rebuild the expression into a power series element."""
        prec = None

        if expr.has(Order):
            var, p = expr.getO().expr.args

            if var != self.symbol:
                raise ValueError("Order contains different symbol than ring.")

            if int(p) > self.prec:
                prec = self.prec
            else:
                prec = int(p)

            expr = expr.removeO()

        def _rebuild(expr):
            if expr.is_Symbol:
                if expr == self.symbol:
                    return self.gen
                else:
                    raise ValueError("Expr contains different symbol than ring.")
            elif expr.is_Add:
                return reduce(add, map(_rebuild, expr.args))
            elif expr.is_Mul:
                return reduce(mul, map(_rebuild, expr.args))
            elif expr.is_Pow:
                base, exp = expr.as_base_exp()
                if exp.is_Integer:
                    try:
                        return _rebuild(base) ** exp
                    except ValueError:
                        raise ValueError("Invalid exponent in power series")
                else:
                    raise ValueError("Invalid exponent in power series")
            else:
                return self.from_ground(self.domain.convert(expr))

        series = _rebuild(sympify(expr))

        if prec is None:
            return series
        else:
            R = self.ring
            s = R(self.to_list(series), prec)
            return self.from_element(s)

    def from_list(
        self, lst: list[Er], prec: int | None = None
    ) -> PowerSeriesElement[Er]:
        """Create a power series element from a list of coefficients."""
        R = self.ring
        s = R.from_list(lst, prec)
        return self.from_element(s)

    def from_element(self, element: TSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Convert a lower power series element to a PowerSeriesElement."""
        return PowerSeriesElement(self, element)

    def from_int(self, arg: int) -> PowerSeriesElement[Er]:
        """Convert an integer to a power series element."""
        g = self.domain_new(arg)
        return self.from_ground(g)

    def from_ground(self, arg: Er) -> PowerSeriesElement[Er]:
        """Convert a ground element to a power series element."""
        R = self.ring
        s = R.from_list([arg])
        return self.from_element(s)

    def to_expr(self, element: PowerSeriesElement[Er]) -> Expr:
        """Convert a power series element to an expression."""
        dom: Domain[Er] = self.domain
        coeffs: list[Er] = self.ring.to_list(element.series)
        prec: int | None = self.ring.series_prec(element.series)
        result: list[Expr] = []
        x: Expr = self.symbol

        if coeffs and coeffs[0] != 0:
            result.append(dom.to_sympy(coeffs[0]))

        for i, coeff in enumerate(coeffs[1:], start=1):
            if coeff == 0:
                continue
            result.append(Mul(dom.to_sympy(coeff), Pow(x, i)))

        if prec is not None:
            result.append(Order(x**prec))

        return Add(*result)

    def to_list(self, element: PowerSeriesElement[Er]) -> list[Er]:
        """Returns a list of coefficients of a power series."""
        return self.ring.to_list(element.series)

    def to_dense(self, element: PowerSeriesElement[Er]) -> dup[Er]:
        """Returns a dense list coefficients of a power series."""
        return self.ring.to_dense(element.series)

    def domain_new(self, arg: Er | int) -> Er:
        """Convert arg to the element of ground domain of ring."""
        return self.domain.convert(arg, self.domain)

    def ring_new(self, arg: Expr | Er | int) -> PowerSeriesElement[Er]:
        """Create a power series element from various types."""
        if isinstance(arg, Expr):
            return self.from_expr(arg)
        elif isinstance(arg, int):
            return self.from_int(arg)
        elif self.domain.of_type(arg):
            return self.from_ground(arg)
        else:
            raise TypeError(f"{arg.__class__.__name__} type not supported")

    __call__ = ring_new

    def square(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """
        Compute the square of a power series.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x", 5)
        >>> x = R.gen
        >>> R.square(1 + x + x**2)
        1 + 2*x + 3*x**2 + 2*x**3 + x**4

        """
        R = self.ring
        series = R.square(s.series)
        return self.from_element(series)

    def compose(
        self, s: PowerSeriesElement[Er], t: PowerSeriesElement[Er]
    ) -> PowerSeriesElement[Er]:
        """
        Compute the composition of two power series.

        The second series must have zero constant term for valid composition.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x", 5)
        >>> x = R.gen

        >>> R.compose(3 + x + 14*x**2, x**2 + x**3)
        3 + x**2 + x**3 + 14*x**4 + O(x**5)

        """
        R = self.ring
        series = R.compose(s.series, t.series)
        return self.from_element(series)

    def inverse(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """
        Compute the multiplicative inverse of a power series.

        The constant term must be a unit in the ground domain for
        valid multiplicative inverse.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x", 5)
        >>> x = R.gen
        >>> R.inverse(1 + x)
        1 - x + x**2 - x**3 + x**4 + O(x**5)

        """
        R = self.ring
        series = R.inverse(s.series)
        return self.from_element(series)

    def reversion(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """
        Compute the compositional inverse of a power series.

        The constant term must be zero and the linear term must be a unit in the ground
        domain for valid compositional inverse.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x", 5)
        >>> x = R.gen
        >>> R.reversion(x + 4*x**2 + 8*x**3)
        x - 4*x**2 + 24*x**3 - 160*x**4 + O(x**5)

        """
        R = self.ring
        series = R.reversion(s.series)
        return self.from_element(series)

    def differentiate(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.series import PowerSeriesRingRing
        >>> R = PowerSeriesRingRing(ZZ, "x", 5)
        >>> x = R.gen
        >>> R.differentiate(10 + x + x**2 + x**3 + x**4 + R.order_term())
        1 + 2*x + 3*x**2 + 4*x**3 + O(x**4)

        """

        R = self.ring
        series = R.differentiate(s.series)
        return self.from_element(series)


class PowerSeriesRingField(PowerSeriesRingRing[Ef], Generic[Ef]):
    """A class for representing Univariate Power Series Rings over a Field."""

    ring: PowerSeriesRingFieldProto[TSeriesElement[Ef], Ef]
    dtype: type[PowerSeriesElement[Ef]]
    domain: Domain[Ef]
    prec: int
    symbol: Expr

    def __init__(self, domain: Domain[Ef], symbol: str | Expr = "x", prec: int = 6):
        if not isinstance(symbol, (str, Expr)):
            raise GeneratorsError("Symbol must be a string or a SymPy expression.")

        if isinstance(symbol, str):
            symbol = Symbol(symbol)

        # The type checker thinks `ring` overrides with an incompatible type
        # because it cannot infer that this attribute is effectively read-only.
        # Safe to ignore here in the future, stricter immutability checks may resolve this.
        self.ring = _power_series_ring(domain, prec)  # type: ignore
        self.dtype = PowerSeriesElement
        self.domain = domain
        self.symbol = symbol
        self.prec = prec

        self.one = self.from_ground(self.domain.one)
        self.zero = self.from_ground(self.domain.zero)
        self.gen = self.from_element(self.ring.gen)

    def sqrt(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the square root of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.sqrt(1 + 2*x + 4*x**2 + 8*x**3)
        1 + x + 3/2*x**2 + 5/2*x**3 - 29/8*x**4 - 1/8*x**5 + O(x**6)

        """
        R = self.ring
        series = R.sqrt(s.series)
        return self.from_element(series)

    def integrate(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the integral of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.integrate(1 + x + x**2 + x**3 + x**4 + x**5 + R.order_term())
        x + 1/2*x**2 + 1/3*x**3 + 1/4*x**4 + 1/5*x**5 + O(x**6)

        """
        R = self.ring
        series = R.integrate(s.series)
        return self.from_element(series)

    def log(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the logarithm of a power series.

        The constant term should be one for proper expansion in the Rational Field.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.log(1 + x)
        x - 1/2*x**2 + 1/3*x**3 - 1/4*x**4 + 1/5*x**5 + O(x**6)

        """
        R = self.ring
        series = R.log(s.series)
        return self.from_element(series)

    def log1p(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the logarithm of a power series plus one.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.log1p(x)
        x - 1/2*x**2 + 1/3*x**3 - 1/4*x**4 + 1/5*x**5 + O(x**6)

        """

        R = self.ring
        series = R.log1p(s.series)
        return self.from_element(series)

    def exp(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the exponential of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.exp(x)
        1 + x + 1/2*x**2 + 1/6*x**3 + 1/24*x**4 + 1/120*x**5 + O(x**6)

        """
        R = self.ring
        series = R.exp(s.series)
        return self.from_element(series)

    def expm1(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the exponential of a power series minus one.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.expm1(x)
        x + 1/2*x**2 + 1/6*x**3 + 1/24*x**4 + 1/120*x**5 + O(x**6)

        """
        R = self.ring
        series = R.expm1(s.series)
        return self.from_element(series)

    def atan(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the arctangent of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.atan(x)
        x - 1/3*x**3 + 1/5*x**5 + O(x**6)

        """
        R = self.ring
        series = R.atan(s.series)
        return self.from_element(series)

    def atanh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the hyperbolic arctangent of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.atanh(x)
        x + 1/3*x**3 + 1/5*x**5 + O(x**6)

        """
        R = self.ring
        series = R.atanh(s.series)
        return self.from_element(series)

    def asin(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the arcsine of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.asin(x)
        x + 1/6*x**3 + 3/40*x**5 + O(x**6)

        """
        R = self.ring
        series = R.asin(s.series)
        return self.from_element(series)

    def asinh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the hyperbolic arcsine of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.asinh(x)
        x - 1/6*x**3 + 3/40*x**5 + O(x**6)

        """
        R = self.ring
        series = R.asinh(s.series)
        return self.from_element(series)

    def tan(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the tangent of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.tan(x)
        x + 1/3*x**3 + 2/15*x**5 + O(x**6)

        """
        R = self.ring
        series = R.tan(s.series)
        return self.from_element(series)

    def tanh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the hyperbolic tangent of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.tanh(x)
        x - 1/3*x**3 + 2/15*x**5 + O(x**6)

        """
        R = self.ring
        series = R.tanh(s.series)
        return self.from_element(series)

    def sin(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the sine of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.sin(x)
        x - 1/6*x**3 + 1/120*x**5 + O(x**6)

        """
        R = self.ring
        series = R.sin(s.series)
        return self.from_element(series)

    def sinh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the hyperbolic sine of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.sinh(x)
        x + 1/6*x**3 + 1/120*x**5 + O(x**6)

        """
        R = self.ring
        series = R.sinh(s.series)
        return self.from_element(series)

    def cos(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the cosine of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.cos(x)
        1 - 1/2*x**2 + 1/24*x**4 + O(x**6)

        """
        R = self.ring
        series = R.cos(s.series)
        return self.from_element(series)

    def cosh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """
        Compute the hyperbolic cosine of a power series.

        The constant term should be zero for proper expansion in the Rational Field.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.cosh(x)
        1 + 1/2*x**2 + 1/24*x**4 + O(x**6)

        """
        R = self.ring
        series = R.cosh(s.series)
        return self.from_element(series)


class PowerSeriesElement(DomainElement, CantSympify, Generic[Er]):
    """A class for representing elements of a Power Series."""

    def __init__(self, ring: PowerSeriesRingRing[Er], series: TSeriesElement[Er]):
        self.series = series
        self._parent_ring = ring
        self.ring = ring.ring

    def __repr__(self) -> str:
        R = self.ring
        return R.pretty(self.series, symbol=str(self._parent_ring.symbol))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PowerSeriesElement):
            return False

        R = self.ring
        return R.equal_repr(self.series, other.series)

    def __hash__(self) -> int:
        return hash((self._parent_ring, self.ring, self.series))

    def __pos__(self) -> PowerSeriesElement[Er]:
        R = self.ring
        s = R.positive(self.series)
        return self._new(s)

    def __neg__(self) -> PowerSeriesElement[Er]:
        R = self.ring
        s = R.negative(self.series)
        return self._new(s)

    def _new(self, series: TSeriesElement[Er]) -> PowerSeriesElement[Er]:
        return self.__class__(self._parent_ring, series)

    def __add__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        """
        Add two power Series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.sin(x) + R.cos(x)
        1 + x - 1/2*x**2 - 1/6*x**3 + 1/24*x**4 + 1/120*x**5 + O(x**6)

        """
        if isinstance(other, PowerSeriesElement):
            if self._parent_ring == other._parent_ring:
                return self._add(other)
            raise ValueError("Cannot add power series with different rings.")

        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._add_ground(domain.convert(other, domain))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            return NotImplemented

    def __radd__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._add_ground(self._parent_ring.domain_new(other))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            return NotImplemented

    def __sub__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        """
        Subtract two power Series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.sin(x) - R.cos(x)
        -1 + x + 1/2*x**2 - 1/6*x**3 - 1/24*x**4 + 1/120*x**5 + O(x**6)

        """
        if isinstance(other, PowerSeriesElement):
            if self._parent_ring == other._parent_ring:
                return self._sub(other)
            raise ValueError("Cannot subtract power series with different rings.")

        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._sub_ground(self._parent_ring.domain_new(other))
        elif domain.of_type(other):
            return self._sub_ground(other)
        else:
            return NotImplemented

    def __rsub__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._rsub_ground(self._parent_ring.domain_new(other))
        elif domain.of_type(other):
            return self._rsub_ground(other)
        else:
            return NotImplemented

    def __mul__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        """
        Multiply two power Series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> R.sin(x) * R.cos(x)
        x - 2/3*x**3 + 2/15*x**5 + O(x**6)

        """
        if isinstance(other, PowerSeriesElement):
            if self._parent_ring == other._parent_ring:
                return self._mul(other)
            raise ValueError("Cannot multiply power series with different rings.")

        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._mul_ground(self._parent_ring.domain_new(other))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            return NotImplemented

    def __rmul__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._mul_ground(domain.convert(other, domain))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            return NotImplemented

    def __pow__(self, n: int) -> PowerSeriesElement[Er]:
        """
        Raise a power Series to a integer power.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> (1 + x)**3
        1 + 3*x + 3*x**2 + x**3

        """
        if n == 0:
            return self._parent_ring.one
        if n == 1:
            return self
        else:
            return self._pow_int(n)

    def __truediv__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        """
        Divide two power Series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series import PowerSeriesRingField
        >>> R = PowerSeriesRingField(QQ, "x")
        >>> x = R.gen
        >>> t = R.sin(x) / R.cos(x); t
        x + 1/3*x**3 + 2/15*x**5 + O(x**6)
        >>> R.tan(x) == t
        True
        """
        if isinstance(other, PowerSeriesElement):
            if self._parent_ring == other._parent_ring:
                return self._div(other)
            raise ValueError("Cannot add power series with different rings.")

        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._div_ground(domain.convert(other, domain))
        elif domain.of_type(other):
            return self._div_ground(other)
        else:
            return NotImplemented

    def __rtruediv__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self._parent_ring.domain
        if isinstance(other, int):
            return self._rdiv_ground(domain.convert(other, domain))
        elif domain.of_type(other):
            return self._rdiv_ground(other)
        else:
            return NotImplemented

    def _add(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.add(self.series, other.series)
        return self._new(series)

    def _add_ground(self, other: Er) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.add_ground(self.series, other)
        return self._new(series)

    def _sub(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.subtract(self.series, other.series)
        return self._new(series)

    def _sub_ground(self, other: Er) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.subtract_ground(self.series, other)
        return self._new(series)

    def _rsub_ground(self, other: Er) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.rsubtract_ground(self.series, other)
        return self._new(series)

    def _mul(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.multiply(self.series, other.series)
        return self._new(series)

    def _mul_ground(self, other: Er) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.multiply_ground(self.series, other)
        return self._new(series)

    def _pow_int(self, n: int) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.pow_int(self.series, n)
        return self._new(series)

    def _div(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring
        series = R.divide(self.series, other.series)
        return self._new(series)

    def _div_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self._parent_ring.from_ground(other)
        return self._div(s2)

    def _rdiv_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self._parent_ring.from_ground(other)
        return s2._div(self)

    def as_expr(self) -> Expr:
        return self._parent_ring.to_expr(self)

    def constant_coefficient(self) -> Er:
        """Return the constant coefficient of the series."""
        R = self.ring
        return R.constant_coefficient(self.series)

    def removeO(self) -> PowerSeriesElement[Er]:
        """Remove the big O notation from the series."""
        R = self.ring
        coeffs = R.to_list(self.series)
        series = R.from_list(coeffs)
        return self._new(series)

    @property
    def is_ground(self) -> bool | None:
        """Check if the series is a ground element."""
        R = self.ring
        return R.is_ground(self.series)
