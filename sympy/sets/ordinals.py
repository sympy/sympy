from sympy.core import Basic
import operator

class Ordinal(Basic):
    """
    Represents ordinals in Cantor normal form.

    Internally, this class is just a list of lists (known
    as "terms").
    In single terms [a, b], a represents power to the ordinal and b represents multiplying factor
    Examples
    ========
    >>> from sympy.sets import Ordinal, Ordinals
    >>> Ordinals.w  # doctest: +SKIP
    w
    >>> a = Ordinals.w
    >>> a.is_limit_ordinal
    True
    >>> Ordinal([[3,2]])  # doctest: +SKIP
    {w**3}*2
    >>> Ordinal([[5,1],[3,2]])  # doctest: +SKIP
    w**5 + {w**3}*2
    """
    def __new__(cls, terms):
        obj = super(Ordinal, cls).__new__(cls, *terms)
        powers = [i[0] for i in obj.args]
        if not all(earlier >= later for earlier, later in zip(powers, powers[1:])):
            raise ValueError("powers must be in decreasing order")
        return obj

    @property
    def is_successor_ordinal(self):
        return (self.args[-1][0] == 0)

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

    def _compare(self, other, index):
        value_self = [i[index] for i in self.args]
        value_other = [i[index] for i in other.args]
        if value_other == value_self:
            return -1
        else:
            return map(operator.eq, value_self, value_other).index(False)

    def __eq__(self, other):
        if isinstance(other, Ordinal):
            return self.args == other.args
        else:
            return False

    def __lt__(self, other):
        if isinstance(other, int):
            return False
        elif isinstance(other, Ordinal):
            index = self._compare(other, 0)
            if index == -1:
                index_factor = self._compare(other, 1)
                if index_factor == -1:
                    return false
                else:
                    return operator.lt(self.args[index_factor][1], other.args[index_factor][1])
            else:
                if index in (len(self.args), len(other.args)):
                    if self.args[index-1][1] == other.args[index-1][1]:
                        return len(self.args) < len(other.args)
                    else:
                        return operator.lt(self.args[index-1][1], other.args[index-1][1])
                else:
                    return operator.lt(self.args[index][0], other.args[index][0])
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
            if i[1] == 0:
                continue
            elif i[0] == 0:
                net_str += str(i[1])
            elif i[0] == 1:
                if i[1] == 1:
                    net_str += 'w'
                else:
                    net_str += 'w*%s'%i[1]
            else :
                if i[1]==1:
                    net_str += 'w**%s'%i[0]
                else:
                    net_str += '{w**%s'%i[0]+'}*%s'%i[1]
            plus_count += 1
        return(net_str)

    __repr__ = __str__

    def __add__(self, other):

        if isinstance(other, int):
            new_terms = list(self.args)
            if other < 0:
                raise ValueError("can only add ordinals and positive integers")

            elif new_terms[-1][0] == 0:
                new_terms[-1][1] += other
                return Ordinal(new_terms)
            else:
                new_terms.append([0, other])
                return Ordinal(new_terms)
        else:
            a_terms = list(self.args)
            b_terms = list(other.args)
            power_self = [i[0] for i in a_terms]
            power_other = [i[0] for i in b_terms]
            b1 = power_other[0]
            r = self._insert_index(power_self, power_other[0])
            if not b1 in power_self:
                a_terms.insert(r,[b1, 0])
            else:
                r-=1
            net = a_terms[:r]
            term2 = [b_terms[0][0], b_terms[0][1]+a_terms[r][1]]
            if term2:
                net.append(term2)
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
        SingletonOrdinal.__instance.w = Ordinal([[1, 1]])
        return SingletonOrdinal.__instance

Ordinals = SingletonOrdinal()
