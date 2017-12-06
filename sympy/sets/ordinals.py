from sympy.core import Basic

def reverse_bisect_right(a, x, lo=0, hi=None):
    """Return the index where to insert item x in list a, assuming a is sorted in descending order

    Essentially, the function returns number of elements in a which are >= than x.
    >>> a = [8, 6, 5, 4, 2]
    >>> reverse_bisect_right(a, 5)
    3
    >>> a[:reverse_bisect_right(a, 5)]
    [8, 6, 5]
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

class Ordinal(Basic):
    """
    Represents ordinals in Cantor normal form.

    Internally, this class is just a list of lists (known
    as "terms").
    In single terms [a,b,c], a represents index of ordinal number
    b represents power to the ordinal and c represents multiplying factor
    Examples
    ========
    >>> from sympy.sets import Ordinal, Ordinals
    >>> Ordinals.w
    w
    >>> Ordinals.w1
    w_1
    >>> a = Ordinals.w
    >>> a.is_limit_ordinal
    True
    >>> Ordinal([[0,3,2]])
    w**3.2
    >>> Ordinal([[1,5,1]])
    w_1**5
    >>> Ordinal([[0,5,1],[0,3,2]])
    w**5 + w**3.2
    """
    def __new__(cls, terms):

        obj = super(Ordinal,cls).__new__(cls,*terms)
        obj.terms = terms
        obj.index = obj.args[0][0]
        return obj

    @property
    def is_countable(self):
        return self.index == 0

    @property
    def is_uncountable(self):
        return self.index > 0

    @property
    def is_successor_ordinal(self):
        return (self.args[-1][1] == 0)

    @property
    def is_limit_ordinal(self):
        return not self.is_successor_ordinal

    def __str__(self):
        def represent(term):
            if term[0] == 0:
                return 'w'
            else:
                return 'w_%s' % term[0]
        net_str = ""
        plus_count=0
        for i in self.args:
            if plus_count:
                net_str+=" + "
            if i[2] == 0:
                continue
            elif i[1] == 0:
                net_str+=str(i[2])
            elif i[1] == 1:
                if i[2] == 1:
                    net_str+=represent(i)
                else:
                    net_str+=represent(i)+'.%s'%i[2]
            else :
                if i[2]==1:
                    net_str+=represent(i)+'**%s'%i[1]
                else:
                    net_str+=represent(i)+'**%s'%i[1]+'.%s'%i[2]
            plus_count += 1
        return(net_str)
    def __repr__(self):
        return str(self)

    def __add__(self,other):

        if isinstance(other, int):
            new_terms = list(self.args)
            if other < 0:
                raise ValueError("can only add ordinals and positive integers")

            elif new_terms[0][1] == 0:
                new_terms[0][2]+=other
                return Ordinal(new_terms)
            else:
                new_terms.append([(self.index), 0, other])
                return Ordinal(new_terms)
        else:
            a_terms = list(self.args)
            b_terms = list(other.args)
            print(a_terms)
            print(b_terms)
            power_self = [i[1] for i in a_terms]
            power_other = [i[1] for i in b_terms]
            b1 = power_other[0]
            r = reverse_bisect_right(power_self, power_other[0])
            print(r)
            if not b1 in power_self:
                a_terms.insert(r,[self.index, b1, 0])
            else:
                r-=1
            print(a_terms)
            net = a_terms[:r]
            print(net)
            term2 = [b_terms[0][0],b_terms[0][1], b_terms[0][2]+a_terms[r][2]]
            print(term2)
            if term2:
                net.append(term2)
            term3 = (b_terms[1:])
            print(term3)
            if term3:
                net.append(*term3)
            print(net)
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
        SingletonOrdinal.__instance.w = Ordinal([[0, 1, 1]])
        # first uncountable ordinal
        SingletonOrdinal.__instance.w1 = Ordinal([[1, 1, 1]])
        return SingletonOrdinal.__instance

Ordinals = SingletonOrdinal()
