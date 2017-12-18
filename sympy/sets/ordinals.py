from sympy.core import Basic
import operator

class OrdinalPow(Basic):
    """
    Represents ordinal exponential and multiplication terms one of the
    building blocks of the Ordinal class.
    In OrdinalPow(a, b) a represents exponent and b represents multiplicity.
    """
    def __new__(cls, *terms):
        if not all (isinstance(i, (int, Ordinal)) for i in terms):
            raise TypeError("arguments must be an integer or an Ordinal instance")
        if not len(terms) == 2:
            raise TypeError("number of arguments is invalid")
        return Basic.__new__(cls, *terms)

    @property
    def exp(self):
        return self.args[0]

    @property
    def mult(self):
        return self.args[1]

    @property
    def contains_ordinal(self):
        return any(isinstance(x, Ordinal) for x in self.args)

    def _compare_term(self, other, op):
        if self.exp == other.exp:
            return op(self.mult, other.mult)
        else:
            return op(self.exp, other.exp)



class Ordinal(Basic):
    """
    Represents ordinals in Cantor normal form.

    Internally, this class is just a list of instances of OrdinalPow
    Examples
    ========
    >>> from sympy.sets import Ordinal, Ordinals, OrdinalPow
    >>> Ordinals.w  # doctest: +SKIP
    w
    >>> a = Ordinals.w
    >>> a.is_limit_ordinal
    True
    >>> Ordinal([OrdinalPow(3,2)])  # doctest: +SKIP
    {w**3}*2
    >>> Ordinal([OrdinalPow(5,1),OrdinalPow(3,2)])  # doctest: +SKIP
    w**5 + {w**3}*2
    """
    def __new__(cls, terms):
        obj = super(Ordinal, cls).__new__(cls, *terms)
        powers = [i.exp for i in obj.args]
        if not all(earlier >= later for earlier, later in zip(powers, powers[1:])):
            raise ValueError("powers must be in decreasing order")
        return obj

    @property
    def is_successor_ordinal(self):
        return (self.args[-1].exp == 0)

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

    def _convert_integer(self, integer_value):
        return Ordinal([OrdinalPow(0, integer_value)])

    def __lt__(self, other):
        if self == other:
            return False

        elif isinstance(other, int):
            return operator.lt(self, self._convert_integer(other))

        elif isinstance(other, Ordinal):
            index =  map(operator.eq, self.args, other.args).index(False)
            if index not in (len(self.args), len(other.args)):
                return self.args[index]._compare_term(other.args[index], operator.lt)
            else:
                return len(self.args) < len(other.args)
        else:
            raise TypeError("cannot compare types: %s and %s" % (type(self), type(other)))

    def __le__(self, other):
        return (self == other or self < other)

    def __gt__(self, other):
        if type(other) not in (int , Ordinal):
            raise TypeError("cannot compare types: %s and %s" % (type(self), type(other)))
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
            elif i.exp == 0:
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

        if isinstance(other, int):
            new_terms = list(self.args)
            if other < 0:
                raise ValueError("can only add ordinals and positive integers")

            elif new_terms[-1].exp == 0:
                new_terms[-1] = OrdinalPow(0, other+new_terms[-1].mult)
                return Ordinal(new_terms)
            else:
                new_terms.append(OrdinalPow(0, other))
                return Ordinal(new_terms)
        else:
            a_terms = list(self.args)
            b_terms = list(other.args)
            power_self = [i.exp for i in a_terms]
            power_other = [i.exp for i in b_terms]
            b1 = power_other[0]
            r = self._insert_index(power_self, power_other[0])
            if not b1 in power_self:
                a_terms.insert(r,OrdinalPow(b1, 0))
            else:
                r-=1
            net = a_terms[:r]
            term2 = [b_terms[0].exp, b_terms[0].mult+a_terms[r].mult]
            if term2:
                net.append(OrdinalPow(*term2))
            term3 = (b_terms[1:])
            if term3:
                net.append(*term3)
            return(Ordinal(net))

    def __radd__(self, other):
        if type(other) is int and other > 0:
            return self
        else:
            raise ValueError("can only add ordinals and positive integers")

class SingletonOrdinal(object):
    __instance = None

    def __new__(cls):
        if SingletonOrdinal.__instance is None:
            SingletonOrdinal.__instance = object.__new__(cls)
        # first countably infinite ordinal
        SingletonOrdinal.__instance.w = Ordinal([OrdinalPow(1, 1)])
        return SingletonOrdinal.__instance

Ordinals = SingletonOrdinal()
