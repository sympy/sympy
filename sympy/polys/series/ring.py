from __future__ import annotations

from typing import Generic, Type, Union
from typing_extensions import TypeIs
from sympy.polys.domains.domain import Domain, Er

from sympy.core.expr import Expr
from sympy.core.symbol import Symbol
from sympy.core.sympify import CantSympify, sympify
from sympy.core.exprtools import decompose_power
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.series.order import Order
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.densebasic import dup_strip, dup
from sympy.polys.series.ringflint import QQSeries, ZZSeries
from sympy.polys.series.ringpython import (
    USeries,
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)
from sympy.external.gmpy import GROUND_TYPES

flint: bool = False

if GROUND_TYPES == "flint":
    flint = True


TSeries = Union[USeries, ZZSeries, QQSeries]


def power_series_ring(
    symbol: str, K: Domain[Er], prec: int = 6
) -> tuple[PowerSeriesRing[Er], PowerSeriesElement[Er]]:
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
    >>> R.sin(x+x**2)
    x + 1/6*x**3 + O(x**4)
    >>> R.cos(x+x**2)
    1 - 1/2*x**2 + O(x**4)
    """

    ring = PowerSeriesRing(K, symbol, prec)
    return ring, ring.gen


class PowerSeriesRing(Generic[Er]):
    """A class for representing Univariate Power Series Rings."""

    is_PowerSeries: bool = True
    has_assoc_Ring: bool = True
    has_assoc_Field: bool = False

    dtype: Type[PowerSeriesElement[Er]]
    domain: Domain[Er]
    symbol: Expr
    prec: int

    def __init__(self, domain: Domain[Er], symbol: str | Expr = "x", prec: int = 6):
        if domain.is_ZZ:
            if flint:
                from sympy.polys.series.ringflint import FlintPowerSeriesRingZZ

                ring = FlintPowerSeriesRingZZ(prec)
            else:
                ring = PythonPowerSeriesRingZZ(prec)
        elif domain.is_QQ:
            if flint:
                from sympy.polys.series.ringflint import FlintPowerSeriesRingQQ

                ring = FlintPowerSeriesRingQQ(prec)
            else:
                ring = PythonPowerSeriesRingQQ(prec)
        else:
            raise ValueError(f"Unsupported ground domain: {domain}")

        if isinstance(symbol, str):
            symbol = Symbol(symbol)

        self.ring = ring
        self.dtype = PowerSeriesElement
        self.domain = domain
        self.symbol = symbol
        self.prec = prec

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

    def is_element(self, element) -> TypeIs[PowerSeriesElement[Er]]:
        """Check if element belongs to this ring."""
        return isinstance(element, PowerSeriesElement) and element.ring == self

    def from_expr(self, expr) -> PowerSeriesElement[Er]:
        poly = self._rebuild_expr(sympify(expr))
        return self.ring_new(poly)

    def _rebuild_expr(self, expr) -> TSeries:
        prec = None
        fring_prec: bool = False

        if expr.has(Order):
            prec = int(expr.getO().expr.args[1])
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

            if pow < len(coeffs):
                coeffs[-pow - 1] = coeff
            else:
                fring_prec = True

        if fring_prec:
            prec = self.prec
        return self.ring(dup_strip(coeffs), prec)

    def from_list(
        self, lst: list[Er], prec: int | None = None
    ) -> PowerSeriesElement[Er]:
        R = self.ring
        s = R.from_list(lst, prec)
        return self.from_element(s)

    def from_element(self, element: TSeries) -> PowerSeriesElement[Er]:
        return PowerSeriesElement(self, element)

    def domain_new(self, arg: Expr | Er | int) -> Er:
        return self.domain.convert(arg, self.domain)

    def ground_new(self, arg: Expr | Er | int) -> PowerSeriesElement[Er]:
        R = self.ring
        g: Er = self.domain_new(arg)
        series = R.from_list([g])
        return self.from_element(series)

    def ring_new(self, arg: TSeries | Expr | Er | int) -> PowerSeriesElement[Er]:
        if isinstance(arg, Expr):
            return self.from_expr(arg)
        elif isinstance(arg, int):
            return self.ground_new(arg)
        elif self.domain.of_type(arg):
            return self.ground_new(arg)
        else:
            return self.from_element(arg)

    __call__ = ring_new

    def square(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the square of a power series."""
        R = self.ring
        series = R.square(s.series)
        return self.from_element(series)

    def sqrt(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the square root of a power series."""
        R = self.ring

        if hasattr(R, "sqrt"):
            series = R.sqrt(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement sqrt."
            )

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

    def integrate(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the integral of a power series."""
        R = self.ring
        if hasattr(R, "integrate"):
            series = R.integrate(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement integrate."
            )

    def log(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the logarithm of a power series."""
        R = self.ring
        if hasattr(R, "log"):
            series = R.log(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement log."
            )

    def logp1(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the logarithm of a power series plus one."""
        R = self.ring
        if hasattr(R, "logp1"):
            series = R.logp1(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement logp1."
            )

    def exp(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the exponential of a power series."""
        R = self.ring
        if hasattr(R, "exp"):
            series = R.exp(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement exp."
            )

    def expm1(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the exponential of a power series minus one."""
        R = self.ring
        if hasattr(R, "expm1"):
            series = R.expm1(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement expm1."
            )

    def atan(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the arctangent of a power series."""
        R = self.ring
        if hasattr(R, "atan"):
            series = R.atan(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement atan."
            )

    def atanh(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the hyperbolic arctangent of a power series."""
        R = self.ring
        if hasattr(R, "atanh"):
            series = R.atanh(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement atanh."
            )

    def asin(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the arcsine of a power series."""
        R = self.ring
        if hasattr(R, "asin"):
            series = R.asin(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement asin."
            )

    def asinh(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the hyperbolic arcsine of a power series."""
        R = self.ring
        if hasattr(R, "asinh"):
            series = R.asinh(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement asinh."
            )

    def tan(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the tangent of a power series."""
        R = self.ring
        if hasattr(R, "tan"):
            series = R.tan(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement tan."
            )

    def tanh(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the hyperbolic tangent of a power series."""
        R = self.ring
        if hasattr(R, "tanh"):
            series = R.tanh(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement tanh."
            )

    def sin(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the sine of a power series."""
        R = self.ring
        if hasattr(R, "sin"):
            series = R.sin(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement sin."
            )

    def sinh(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the hyperbolic sine of a power series."""
        R = self.ring
        if hasattr(R, "sinh"):
            series = R.sinh(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement sinh."
            )

    def cos(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the cosine of a power series."""
        R = self.ring
        if hasattr(R, "cos"):
            series = R.cos(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement cos."
            )

    def cosh(self, s: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        """Return the hyperbolic cosine of a power series."""
        R = self.ring
        if hasattr(R, "cosh"):
            series = R.cosh(s.series)
            return self.from_element(series)
        else:
            raise NotImplementedError(
                f"Power Series Ring over ground domain {self.domain} does not implement cosh."
            )


class PowerSeriesElement(DomainElement, CantSympify, Generic[Er]):
    """A class for representing elements of a Power Series."""

    def __init__(self, ring: PowerSeriesRing[Er], series: TSeries):
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

    def new(self, series: TSeries) -> PowerSeriesElement[Er]:
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

    @property
    def is_ground(self) -> bool:
        """Check if the series is a ground element."""
        R = self.ring.ring
        return len(R.to_list(self.series)) <= 1
