from __future__ import annotations
from typing import Any, Callable, Tuple, Union, Optional, cast
from sympy.logic.boolalg import Boolean

from sympy.core import Basic, Integer
import operator


class OmegaPower(Basic):
    """
    Represents ordinal exponential and multiplication terms one of the
    building blocks of the :class:`Ordinal` class.
    In ``OmegaPower(a, b)``, ``a`` represents exponent and ``b`` represents multiplicity.
    """
    def __new__(cls, a: Any, b: Any) -> OmegaPower:
        if isinstance(b, int):
            b = Integer(b)
        if not isinstance(b, Integer) or b <= 0:
            raise TypeError("multiplicity must be a positive integer")

        if not isinstance(a, Ordinal):
            a = Ordinal.convert(a)

        # Cast required: Basic.__new__ returns Basic, but mypy expects OmegaPower
        return cast(OmegaPower, Basic.__new__(cls, a, b))

    @property
    def exp(self) -> Ordinal:
        # Cast required: args[0] is strictly stored as Basic, but logic ensures it's Ordinal
        return cast(Ordinal, self.args[0])

    @property
    def mult(self) -> Integer:
        # Cast required: args[1] is stored as Basic
        return cast(Integer, self.args[1])

    def _compare_term(self, other: OmegaPower, op: Callable[[Any, Any], bool]) -> bool:
        if self.exp == other.exp:
            return op(self.mult, other.mult)
        else:
            return op(self.exp, other.exp)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OmegaPower):
            try:
                other = OmegaPower(0, other)
            except TypeError:
                return NotImplemented
        return self.args == other.args # type: ignore

    def __hash__(self) -> int:
        return Basic.__hash__(self)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, OmegaPower):
            try:
                other = OmegaPower(0, other)
            except TypeError:
                return NotImplemented  # type: ignore
        return self._compare_term(other, operator.lt)


class Ordinal(Basic):
    """
    Represents ordinals in Cantor normal form.

    Internally, this class is just a list of instances of OmegaPower.

    Examples
    ========
    >>> from sympy import Ordinal, OmegaPower
    >>> from sympy.sets.ordinals import omega
    >>> w = omega
    >>> w.is_limit_ordinal
    True
    >>> Ordinal(OmegaPower(w + 1, 1), OmegaPower(3, 2))
    w**(w + 1) + w**3*2
    >>> 3 + w
    w
    >>> (w + 1) * w
    w**2

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Ordinal_arithmetic
    """
    def __new__(cls, *terms: OmegaPower) -> Ordinal:
        obj = super().__new__(cls, *terms)
        powers = [i.exp for i in obj.args]
        if not all(powers[i] >= powers[i+1] for i in range(len(powers) - 1)):
            raise ValueError("powers must be in decreasing order")
        
        # Cast required: super().__new__ returns Basic
        return cast(Ordinal, obj)

    @property
    def terms(self) -> Tuple[OmegaPower, ...]:
        return self.args  # type: ignore

    @property
    def leading_term(self) -> OmegaPower:
        if self == ord0:
            raise ValueError("ordinal zero has no leading term")
        return self.terms[0]

    @property
    def trailing_term(self) -> OmegaPower:
        if self == ord0:
            raise ValueError("ordinal zero has no trailing term")
        return self.terms[-1]

    @property
    def is_successor_ordinal(self) -> bool:
        try:
            return self.trailing_term.exp == ord0
        except ValueError:
            return False

    @property
    def is_limit_ordinal(self) -> bool:
        try:
            return not self.trailing_term.exp == ord0
        except ValueError:
            return False

    @property
    def degree(self) -> Ordinal:
        return self.leading_term.exp

    @classmethod
    def convert(cls, integer_value: Any) -> Ordinal:
        if integer_value == 0:
            return ord0
        return Ordinal(OmegaPower(0, integer_value))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return self.terms == other.terms

    def __hash__(self) -> int:
        return hash(self.args)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented # type: ignore
        
        for term_self, term_other in zip(self.terms, other.terms): # type: ignore
            if term_self != term_other:
                return term_self < term_other
        return len(self.terms) < len(other.terms) # type: ignore

    def __le__(self, other: object) -> bool:
        return (self == other or self < other)

    def __gt__(self, other: object) -> bool:
        return not self <= other

    def __ge__(self, other: object) -> bool:
        return not self < other

    def __str__(self) -> str:
        net_str = ""
        if self == ord0:
            return 'ord0'
        for plus_count, i in enumerate(self.terms):
            if plus_count:
                net_str += " + "

            if i.exp == ord0:
                net_str += str(i.mult)
            elif i.exp == 1:
                net_str += 'w'
            elif isinstance(i.exp, Ordinal) and (len(i.exp.terms) > 1 or i.exp.is_limit_ordinal):
                net_str += 'w**(%s)' % i.exp
            else:
                net_str += 'w**%s' % i.exp

            if not i.mult == 1 and not i.exp == ord0:
                net_str += '*%s' % i.mult
        return net_str

    __repr__ = __str__

    def __add__(self, other: object) -> Ordinal | Any:
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        if other == ord0:
            return self
        a_terms = list(self.terms)
        b_terms = list(other.terms)
        r = len(a_terms) - 1
        b_exp = other.degree
        while r >= 0 and a_terms[r].exp < b_exp:
            r -= 1
        if r < 0:
            terms = b_terms
        elif a_terms[r].exp == b_exp:
            sum_term = OmegaPower(b_exp, a_terms[r].mult + other.leading_term.mult)
            terms = a_terms[:r] + [sum_term] + b_terms[1:]
        else:
            terms = a_terms[:r+1] + b_terms
        return Ordinal(*terms)

    def __radd__(self, other: object) -> Ordinal | Any:
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return other + self

    def __mul__(self, other: object) -> Ordinal | Any:
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        if ord0 in (self, other):
            return ord0
        a_exp = self.degree
        a_mult = self.leading_term.mult
        summation = []
        if other.is_limit_ordinal:
            for arg in other.terms:
                summation.append(OmegaPower(a_exp + arg.exp, arg.mult)) # type: ignore

        else:
            for arg in other.terms[:-1]:
                summation.append(OmegaPower(a_exp + arg.exp, arg.mult)) # type: ignore
            b_mult = other.trailing_term.mult
            summation.append(OmegaPower(a_exp, a_mult*b_mult)) # type: ignore
            summation += list(self.terms[1:])
        return Ordinal(*summation)

    def __rmul__(self, other: object) -> Ordinal | Any:
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return other * self

    def __pow__(self, other: Any) -> Ordinal | Any:
        if not self == omega:
            return NotImplemented
        return Ordinal(OmegaPower(other, 1))


class OrdinalZero(Ordinal):
    """The ordinal zero."""
    pass


class OrdinalOmega(Ordinal):
    """The ordinal omega."""
    def __new__(cls) -> OrdinalOmega:
        # Cast required
        return cast(OrdinalOmega, Ordinal.__new__(cls))

    @property
    def terms(self) -> Tuple[OmegaPower, ...]:
        return (OmegaPower(1, 1),)


ord0 = OrdinalZero()
omega = OrdinalOmega()