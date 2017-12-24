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
        if not isinstance(b, (int, Integer)) or b <= 0:
            raise TypeError("multiplicity must be a positive integer")
        if a == 0:
            a = Ord0
        elif not isinstance(a, Ordinal):
            try:
                a = Ordinal.convert(a)
            except TypeError:
                return NotImplemented

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
    >>> from sympy.sets import Ordinal, Ordinals, OmegaPower
    >>> Ordinals.w  # doctest: +SKIP
    w
    >>> a = Ordinals.w
    >>> a.is_limit_ordinal
    True
    >>> Ordinal(OmegaPower(3,2))  # doctest: +SKIP
    {w**3}*2
    >>> Ordinal(OmegaPower(5,1),OmegaPower(3,2))  # doctest: +SKIP
    w**5 + {w**3}*2
    """
    def __new__(cls, *terms):
        obj = super(Ordinal, cls).__new__(cls, *terms)
        powers = [i.exp for i in obj.args]
        if not all(earlier >= later for earlier, later in zip(powers, powers[1:])):
            raise ValueError("powers must be in decreasing order")
        return obj

    @property
    def is_successor_ordinal(self):
        return self.args[-1].exp == 0

    @property
    def is_limit_ordinal(self):
        return not self.is_successor_ordinal

    @classmethod
    def convert(cls, integer_value):
        return Ordinal(OmegaPower(0, integer_value))

    def __eq__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return self.args == other.args

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
        for i in self.args:
            if plus_count:
                net_str += " + "
            if i.mult == 0:
                continue
            elif i.exp == Ord0:
                net_str += str(i.mult)
            elif i.exp == 1:
                if i.mult == 1:
                    net_str += 'w'
                else:
                    net_str += 'w*%s'%i.mult
            else :
                if i.mult==1:
                    net_str += 'w**%s'%i.exp
                else:
                    net_str += '{w**%s'%i.exp+'}*%s'%i.mult
            plus_count += 1
        return(net_str)



    __repr__ = __str__

    def __add__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        if other == Ord0:
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
        if Ord0 in (self, other):
            return Ord0

        a1 = self.args[0].exp
        prod = []
        for arg in other.args:
            prod.append(OmegaPower(a1 + arg.exp, arg.mult))
        return Ordinal(*prod)

    def __rmul__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return other * self


class Ord0(with_metaclass(Singleton, Ordinal)):

    def __new__(cls):
        return Ordinal.__new__(cls)

class OrdOmega(with_metaclass(Singleton, Ordinal)):

    def __new__(cls):
        return Ordinal.__new__(cls, OmegaPower(1,1))

Ord0 = S.Ord0
OrdOmega = S.OrdOmega
