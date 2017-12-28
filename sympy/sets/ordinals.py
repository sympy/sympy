from sympy.core import Basic, Integer
from sympy.core.compatibility import with_metaclass
from sympy.core.singleton import Singleton, S
import operator

class OmegaPower(Basic):
    """
    Represents ordinal exponential and multiplication terms one of the
    building blocks of the Ordinal class.
    In OmegaPower(a, b) a represents exponent and b represents multiplicity.
    """
    def __new__(cls, a, b):
        if isinstance(b, int):
            b = Integer(b)
        if not isinstance(b, Integer) or b <= 0:
            raise TypeError("multiplicity must be a positive integer")

        if not isinstance(a, Ordinal):
            a = Ordinal.convert(a)

        return Basic.__new__(cls, a, b)

    @property
    def exp(self):
        return self.args[0]

    @property
    def mult(self):
        return self.args[1]

    def _compare_term(self, other, op):
        if self.exp == other.exp:
            return op(self.mult, other.mult)
        else:
            return op(self.exp, other.exp)

    def __eq__(self, other):
        if not isinstance(other, OmegaPower):
            try:
                other = OmegaPower(0, other)
            except TypeError:
                return NotImplemented
        return hash(self) == hash(other)

    def __lt__(self, other):
        if not isinstance(other, OmegaPower):
            try:
                other = OmegaPower(0, other)
            except TypeError:
                return NotImplemented
        return self._compare_term(other, operator.lt)

    def __hash__(self):
        return hash((self.exp, self.mult))

class Ordinal(Basic):
    """
    Represents ordinals in Cantor normal form.

    Internally, this class is just a list of instances of OmegaPower
    Examples
    ========
    >>> from sympy.sets import Ordinal, ord0, OmegaPower
    >>> from sympy.sets.ordinals import omega
    >>> w = omega
    >>> w.is_limit_ordinal
    True
    >>> Ordinal(OmegaPower(w + 1 ,1), OmegaPower(3, 2)) #doctest: +SKIP
    w**(w + 1) + w**3*2
    >>> 3 + w  #doctest: +SKIP
    w
    >>> (w + 1) * w #doctest: +SKIP
    w**2
    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Ordinal_arithmetic
    """
    def __new__(cls, *terms):
        obj = super(Ordinal, cls).__new__(cls, *terms)
        powers = [i.exp for i in obj.args]
        if not all(earlier >= later for earlier, later in zip(powers, powers[1:])):
            raise ValueError("powers must be in decreasing order")
        return obj

    @property
    def is_successor_ordinal(self):
        return self.args[-1].exp == ord0

    @property
    def is_limit_ordinal(self):
        return not self.is_successor_ordinal

    @classmethod
    def convert(cls, integer_value):
        if integer_value == 0:
            return ord0
        return Ordinal(OmegaPower(0, integer_value))

    def __eq__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return self.args == other.args

    def __hash__(self):
        return hash(self.args)

    def __lt__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        for term_self, term_other in zip(self.args, other.args):
            if term_self != term_other:
                return term_self < term_other
        return len(self.args) < len(other.args)

    def __le__(self, other):
        return (self == other or self < other)

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return not self < other

    def __str__(self):
        net_str = ""
        plus_count = 0
        if self == ord0:
            return 'ord0'
        for i in self.args:
            if plus_count:
                net_str += " + "

            if i.exp == ord0:
                net_str += str(i.mult)
            elif i.exp == 1:
                net_str += 'w'
            elif len(i.exp.args) > 1 or i.exp.is_limit_ordinal:
                net_str += 'w**(%s)'%i.exp
            else:
                net_str += 'w**%s'%i.exp

            if not i.mult == 1 and not i.exp == ord0:
                net_str += '*%s'%i.mult

            plus_count += 1
        return(net_str)

    __repr__ = __str__

    def __add__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        if other == ord0:
            return self
        a_terms = list(self.args)
        b_terms = list(other.args)
        r = len(a_terms) - 1
        b_exp = b_terms[0].exp
        while r >= 0 and a_terms[r].exp < b_exp:
            r -= 1
        if r < 0:
            terms = b_terms
        elif a_terms[r].exp == b_exp:
            sum_term = OmegaPower(b_exp, a_terms[r].mult + b_terms[0].mult)
            terms = a_terms[:r] + [sum_term] + b_terms[1:]
        else:
            terms = a_terms[:r+1] + b_terms
        return Ordinal(*terms)

    def __radd__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return other + self

    def __mul__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        if ord0 in (self, other):
            return ord0
        a_exp = self.args[0].exp
        a_mult = self.args[0].mult
        sum = []
        if other.is_limit_ordinal:
            for arg in other.args:
                sum.append(OmegaPower(a_exp + arg.exp, arg.mult))

        else:
            for arg in other.args[:-1]:
                sum.append(OmegaPower(a_exp + arg.exp, arg.mult))
            b_mult = other.args[-1].mult
            sum.append(OmegaPower(a_exp, a_mult*b_mult))
            sum += list(self.args[1:])
        return Ordinal(*sum)

    def __rmul__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return other * self

    def __pow__(self, other):
        if not self == omega:
            return NotImplemented
        return Ordinal(OmegaPower(other, 1))

class OrdinalZero(with_metaclass(Singleton, Ordinal)):
    """The ordinal zero.

    OrdinalZero is a singleton, and can be accessed by ``S.OrdinalZero``
    or can be imported as ``ord0``.

    Examples
    ========

    >>> from sympy import ord0, S
    >>> ord0 is S.OrdinalZero
    True
    """
    pass

class OrdinalOmega(with_metaclass(Singleton, Ordinal)):
    """The ordinal omega which forms the base of all ordinals in cantor normal form.

    OrdinalOmega is a singleton, and can be accessed by ``S.OrdinalOmega``
    or can be imported as ``omega``.

    Examples
    ========

    >>> from sympy.sets.ordinals import omega
    >>> from sympy import S
    >>> omega is S.OrdinalOmega
    True
    """
    def __new__(cls):
        return Ordinal.__new__(cls, OmegaPower(1, 1))

ord0 = S.OrdinalZero
omega = S.OrdinalOmega
