from __future__ import annotations

from sympy.core.expr import Expr
from sympy.core.symbol import Symbol
from sympy.core.exprtools import decompose_power
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.core.sympify import CantSympify, sympify
from sympy.polys.domains.domain import Er, Ef, Domain, DomainElement
from sympy.polys.densebasic import dup
from sympy.polys.series.base import PowerSeriesRingProto, PowerSeriesRingFieldProto
from sympy.polys.series.tring import TElement, _power_series_ring
from sympy.series.order import Order


from typing import Generic, cast, overload


@overload
def power_series_ring(
    symbol: str, K: Domain[Ef], prec: int = 6
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

    is_PowerSeries: bool = True
    has_assoc_Ring: bool = True
    has_assoc_Field: bool = False

    dtype: type[PowerSeriesElement[Er]]
    domain: Domain[Er]
    symbol: Expr
    prec: int

    def __init__(self, domain: Domain[Er], symbol: str | Expr = "x", prec: int = 6):
        if isinstance(symbol, str):
            symbol = Symbol(symbol)

        self.dtype = PowerSeriesElement
        self.domain = domain
        self.symbol = symbol
        self.prec = prec

    @property
    def ring(self) -> PowerSeriesRingProto[TElement, Er]:
        """Return the base ring of the power series."""
        return cast(
            "PowerSeriesRingProto[TElement, Er]",
            _power_series_ring(self.domain, self.prec),
        )

    @property
    def one(self) -> PowerSeriesElement[Er]:
        return self.ground_new(self.domain.one)

    @property
    def zero(self) -> PowerSeriesElement[Er]:
        return self.ground_new(self.domain.zero)

    @property
    def gen(self) -> PowerSeriesElement[Er]:
        return self.from_element(self.ring.gen)

    def __repr__(self) -> str:
        return f"{self.domain}[[{self.symbol}], {self.prec}]"

    def __eq__(self, other) -> bool:
        return self.ring == other.ring

    def is_element(self, element) -> bool:
        """Check if element belongs to this ring."""
        return isinstance(element, PowerSeriesElement) and element.ring == self

    def order_term(self):
        o = self.ring([], self.prec)
        return self.from_element(o)

    def from_expr(self, expr) -> PowerSeriesElement[Er]:
        poly = self._rebuild_expr(sympify(expr))
        return self.ring_new(poly)

    def _rebuild_expr(self, expr) -> TElement:
        prec = None
        fring_prec: bool = False

        if expr.has(Order):
            var, p = expr.getO().expr.args

            if var != self.symbol:
                raise ValueError("Order contains different symbol than ring.")

            if int(p) > self.prec:
                prec = self.prec
            else:
                prec = int(p)

            expr = expr.removeO()

        max_prec = prec if prec else self.prec
        coeffs = [self.domain.zero] * max_prec

        for term in Add.make_args(expr):
            coeff = self.domain.one
            pow = 0

            for factor in Mul.make_args(term):
                if factor.is_Number:
                    coeff = self.domain_new(factor)
                else:
                    base, exp = decompose_power(factor)
                    if base == self.symbol:
                        pow = int(exp)
                    else:
                        raise ValueError("Expr contains different symbol than ring.")

            if pow < max_prec:
                coeffs[pow] = coeff
            else:
                fring_prec = True

        if fring_prec:
            prec = self.prec

        while coeffs and coeffs[-1] == self.domain.zero:
            coeffs.pop()

        return self.ring.from_list(coeffs, prec)

    def from_list(
        self, lst: list[Er], prec: int | None = None
    ) -> PowerSeriesElement[Er]:
        R = self.ring
        s = R.from_list(lst, prec)
        return self.from_element(s)

    def from_element(self, element: TElement) -> PowerSeriesElement[Er]:
        return PowerSeriesElement(self, element)

    def to_list(self, element: TElement | PowerSeriesElement[Er]) -> list[Er]:
        if isinstance(element, PowerSeriesElement):
            return self.ring.to_list(element.series)
        return self.ring.to_list(element)

    def to_dense(self, element: TElement | PowerSeriesElement[Er]) -> dup[Er]:
        if isinstance(element, PowerSeriesElement):
            return self.ring.to_dense(element.series)
        return self.ring.to_dense(element)

    def domain_new(self, arg: Expr | Er | int) -> Er:
        return self.domain.convert(arg, self.domain)

    def ground_new(self, arg: Expr | Er | int) -> PowerSeriesElement[Er]:
        R = self.ring
        g: Er = self.domain_new(arg)
        series = R.from_list([g])
        return self.from_element(series)

    def ring_new(self, arg: TElement | Expr | Er | int) -> PowerSeriesElement[Er]:
        if isinstance(arg, Expr):
            return self.from_expr(arg)
        elif isinstance(arg, int):
            return self.ground_new(arg)
        elif self.ring.is_element(arg):
            return self.from_element(arg)
        elif self.domain.of_type(arg):
            return self.ground_new(arg)
        else:
            raise NotImplementedError

    __call__ = ring_new

    def square(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the square of a power series."""
        R = self.ring
        series = R.square(s.series)
        return self.from_element(series)

    def compose(
        self, s: PowerSeriesElement[Er], t: PowerSeriesElement[Er]
    ) -> PowerSeriesElement[Er]:
        """Return the composition of two power series."""
        R = self.ring
        series = R.compose(s.series, t.series)
        return self.from_element(series)

    def inverse(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the inverse of a power series."""
        R = self.ring
        series = R.inverse(s.series)
        return self.from_element(series)

    def reversion(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the reversion of a power series."""
        R = self.ring
        series = R.reversion(s.series)
        return self.from_element(series)

    def differentiate(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the derivative of a power series."""
        R = self.ring
        series = R.differentiate(s.series)
        return self.from_element(series)


class PowerSeriesRingField(PowerSeriesRingRing[Ef], Generic[Ef]):
    dtype: type[PowerSeriesElement[Ef]]
    domain: Domain[Ef]
    prec: int
    symbol: Expr

    def __init__(self, domain: Domain[Ef], symbol: str | Expr = "x", prec: int = 6):
        if isinstance(symbol, str):
            symbol = Symbol(symbol)

        self.dtype = PowerSeriesElement
        self.domain = domain
        self.symbol = symbol
        self.prec = prec

    @property
    def ring(self) -> PowerSeriesRingFieldProto[TElement, Ef]:
        """Return the base ring of the power series."""
        return cast(
            "PowerSeriesRingFieldProto[TElement, Ef]",
            _power_series_ring(self.domain, self.prec),
        )

    def sqrt(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the square root of a power series."""
        R = self.ring
        series = R.sqrt(s.series)
        return self.from_element(series)

    def integrate(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the integral of a power series."""
        R = self.ring
        series = R.integrate(s.series)
        return self.from_element(series)

    def log(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the logarithm of a power series."""
        R = self.ring
        series = R.log(s.series)
        return self.from_element(series)

    def log1p(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the logarithm of a power series plus one."""
        R = self.ring
        series = R.log1p(s.series)
        return self.from_element(series)

    def exp(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the exponential of a power series."""
        R = self.ring
        series = R.exp(s.series)
        return self.from_element(series)

    def expm1(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the exponential of a power series minus one."""
        R = self.ring
        series = R.expm1(s.series)
        return self.from_element(series)

    def atan(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the arctangent of a power series."""
        R = self.ring
        series = R.atan(s.series)
        return self.from_element(series)

    def atanh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the hyperbolic arctangent of a power series."""
        R = self.ring
        series = R.atanh(s.series)
        return self.from_element(series)

    def asin(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the arcsine of a power series."""
        R = self.ring
        series = R.asin(s.series)
        return self.from_element(series)

    def asinh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the hyperbolic arcsine of a power series."""
        R = self.ring
        series = R.asinh(s.series)
        return self.from_element(series)

    def tan(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the tangent of a power series."""
        R = self.ring
        series = R.tan(s.series)
        return self.from_element(series)

    def tanh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the hyperbolic tangent of a power series."""
        R = self.ring
        series = R.tanh(s.series)
        return self.from_element(series)

    def sin(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the sine of a power series."""
        R = self.ring
        series = R.sin(s.series)
        return self.from_element(series)

    def sinh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the hyperbolic sine of a power series."""
        R = self.ring
        series = R.sinh(s.series)
        return self.from_element(series)

    def cos(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the cosine of a power series."""
        R = self.ring
        series = R.cos(s.series)
        return self.from_element(series)

    def cosh(self, s: PowerSeriesElement[Ef]) -> PowerSeriesElement[Ef]:
        """Return the hyperbolic cosine of a power series."""
        R = self.ring
        series = R.cosh(s.series)
        return self.from_element(series)


class PowerSeriesElement(DomainElement, CantSympify, Generic[Er]):
    """A class for representing elements of a Power Series."""

    def __init__(self, ring: PowerSeriesRingRing[Er], series: TElement):
        self.series = series
        self.ring = ring

    def __repr__(self) -> str:
        R = self.ring.ring
        return R.pretty(self.series, symbol=str(self.ring.symbol))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PowerSeriesElement):
            return False

        R = self.ring.ring
        return R.equal_repr(self.series, other.series)

    def __hash__(self) -> int:
        return hash((self.ring, self.series))

    def __pos__(self) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        s = R.positive(self.series)
        return self.new(s)

    def __neg__(self) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        s = R.negative(self.series)
        return self.new(s)

    def new(self, series: TElement) -> PowerSeriesElement[Er]:
        return self.__class__(self.ring, series)

    def __add__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        if not other:
            return self

        if isinstance(other, PowerSeriesElement):
            return self._add(other)

        domain = self.ring.domain
        if isinstance(other, int):
            return self._add_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            raise NotImplementedError

    def __radd__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._add_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._add_ground(other)
        else:
            raise NotImplementedError

    def __sub__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        if not other:
            return self

        if isinstance(other, PowerSeriesElement):
            return self._sub(other)

        domain = self.ring.domain
        if isinstance(other, int):
            return self._sub_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._sub_ground(other)
        else:
            raise NotImplementedError

    def __rsub__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._rsub_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._rsub_ground(other)
        else:
            raise NotImplementedError

    def __mul__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        if not other:
            return self

        if isinstance(other, PowerSeriesElement):
            return self._mul(other)

        domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            raise NotImplementedError

    def __rmul__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._mul_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._mul_ground(other)
        else:
            raise NotImplementedError

    def __pow__(self, n: int) -> PowerSeriesElement[Er]:
        if n == 0:
            return self.ring.zero
        if n == 1:
            return self
        else:
            return self._pow_int(n)

    def __truediv__(
        self, other: PowerSeriesElement[Er] | Er | int
    ) -> PowerSeriesElement[Er]:
        if not other:
            return self

        if isinstance(other, PowerSeriesElement):
            return self._div(other)

        domain = self.ring.domain
        if isinstance(other, int):
            return self._div_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._div_ground(other)
        else:
            raise NotImplementedError

    def __rtruediv__(self, other: Er | int) -> PowerSeriesElement[Er]:
        domain = self.ring.domain
        if isinstance(other, int):
            return self._rdiv_ground(self.ring.domain_new(other))
        elif domain.of_type(other):
            return self._rdiv_ground(other)
        else:
            raise NotImplementedError

    def _add(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        series = R.add(self.series, other.series)
        return self.new(series)

    def _add_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self.ring.ground_new(other)
        return self._add(s2)

    def _sub(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        series = R.subtract(self.series, other.series)
        return self.new(series)

    def _sub_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self.ring.ground_new(other)
        return self._sub(s2)

    def _rsub_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self.ring.ground_new(other)
        return s2._sub(self)

    def _mul(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        series = R.multiply(self.series, other.series)
        return self.new(series)

    def _mul_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self.ring.ground_new(other)
        return self._mul(s2)

    def _pow_int(self, n: int) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        series = R.pow_int(self.series, n)
        return self.new(series)

    def _div(self, other: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        R = self.ring.ring
        series = R.divide(self.series, other.series)
        return self.new(series)

    def _div_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self.ring.ground_new(other)
        return self._div(s2)

    def _rdiv_ground(self, other: Er) -> PowerSeriesElement[Er]:
        s2 = self.ring.ground_new(other)
        return s2._div(self)

    def to_list(self) -> list[Er]:
        R = self.ring.ring
        return R.to_list(self.series)

    def to_dense(self) -> dup[Er]:
        """Convert the power series to a dense representation."""
        R = self.ring.ring
        return R.to_dense(self.series)

    def as_expr(self) -> Expr:
        dom: Domain[Er] = self.ring.domain
        coeffs: list[Er] = self.ring.ring.to_list(self.series)
        prec: int | None = self.ring.ring.series_prec(self.series)
        result: list[Expr] = []
        x: Expr = self.ring.symbol

        if coeffs and coeffs[0] != 0:
            result.append(dom.to_sympy(coeffs[0]))

        for i, coeff in enumerate(coeffs[1:], start=1):
            if coeff == 0:
                continue
            result.append(Mul(dom.to_sympy(coeff), Pow(x, i)))

        if prec is not None:
            result.append(Order(x**prec))

        return Add(*result)

    def leading_coefficient(self) -> Er:
        """Return the leading coefficient of the series."""
        R = self.ring.ring
        coeffs = R.to_list(self.series)
        if len(coeffs) == 0:
            return self.ring.domain.zero
        return self.ring.domain_new(coeffs[0])

    def removeO(self) -> PowerSeriesElement[Er]:
        """Remove the big O notation from the series."""
        R = self.ring.ring
        coeffs = R.to_list(self.series)
        series = R.from_list(coeffs)
        return self.new(series)

    @property
    def is_ground(self) -> bool:
        """Check if the series is a ground element."""
        R = self.ring.ring
        return len(R.to_list(self.series)) <= 1
