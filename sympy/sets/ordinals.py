from sympy.core import Basic, Integer
import operator

class OmegaPower(Basic):
    """
    Represents ordinal exponential and multiplication terms one of the
    building blocks of the Ordinal class.
    In OmegaPower(a, b) a represents exponent and b represents multiplicity.
    """
    def __new__(cls, a, b):
        if b < 0 or not isinstance(b, (int, Integer)):
            raise TypeError("multiplicity must be a positive integer")
        if a == 0:
            a = Ordinal()
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

    def _insert_index(self, a, x, lo=0, hi=None):
        """Return the index where to insert item x in list a,
        assuming a is sorted in descending order

        Essentially, the function returns number of elements in a which are >= than x.
        """
        if lo < 0:
            raise ValueError('lo must be non-negative')
        if hi is None:
            hi = len(a)
        while lo < hi:
            mid = (lo+hi)//2
            if x > a[mid]: hi = mid
            else: lo = mid+1
        return lo

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
            elif i.exp == Ordinal():
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

        a_terms = list(self.args)
        b_terms = list(other.args)
        power_self = [i.exp for i in a_terms]
        power_other = [i.exp for i in b_terms]
        b1 = power_other[0]
        r = self._insert_index(power_self, power_other[0])
        if not b1 in power_self:
            a_terms.insert(r,OmegaPower(b1, 0))
        else:
            r-=1
        net = a_terms[:r]
        term2 = [b_terms[0].exp, b_terms[0].mult+a_terms[r].mult]
        if term2:
            net.append(OmegaPower(*term2))
        term3 = (b_terms[1:])
        if term3:
            net.append(*term3)
        return(Ordinal(*net))

    def __radd__(self, other):
        if not isinstance(other, Ordinal):
            try:
                other = Ordinal.convert(other)
            except TypeError:
                return NotImplemented
        return other + self

class SingletonOrdinal(object):
    __instance = None

    def __new__(cls):
        if SingletonOrdinal.__instance is None:
            SingletonOrdinal.__instance = object.__new__(cls)
        # first countably infinite ordinal
        SingletonOrdinal.__instance.w = Ordinal(OmegaPower(1, 1))
        return SingletonOrdinal.__instance

Ordinals = SingletonOrdinal()
