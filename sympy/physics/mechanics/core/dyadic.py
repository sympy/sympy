#Note - This is a WIP
#For the time being, Symbols have been used in place of BaseVectors
#to test whether the methods behave as expected.
#Some example usage at - http://pastebin.com/WxJ15ewk

from sympy.core import Expr, Mul, Add, Pow, S, sympify
from sympy.core.decorators import call_highest_priority
#from sympy.vector import BaseVector, VectAdd, VectMul, Vector
#from sympy.physics.mechanics import _check_vector


class Dyadic(Expr):
    """ Dyadic super class """
    _op_priority = 11.0

    def __neg__(self):
        return DyadicMul(S(-1), self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return DyadicAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return DyadicAdd(other, self)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return DyadicAdd(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return DyadicAdd(other, -self)
    
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return DyadicMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return DyadicMul(other, self)

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _dyad_div(self, other)

    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise TypeError("Cannot divide by dyadic")

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def factor(self, *args, **kwargs):
        raise TypeError("Factoring not supported for dyadics")

    def evalf(self, *args):
        return self

    def simplify(self):
        simplify_components = {}
        for x in self.components:
            simplify_components[x] = simplify(self.components[x])
        simplify_components = [x * simplify_components[x] \
                               for x in simplify_components]
        return DyadicAdd(*simplify_components)

    @call_highest_priority('__rand__')
    def __and__(self, other):
        if isinstance(other, Symbol):
            outvec = 0
            for x in self.components:
                vect_dot = Symbol("(" + str(x.args[1]) + "." + str(other) + ")")
                outvec += vect_dot * self.components[x] * x.args[0]
            return outvec
        elif isinstance(other, Dyadic):
            outdyad = 0
            for x in self.components:
                for y in other.components:
                    vect_dot = Symbol("(" + str(x.args[1]) + "." + str(y.args[0]) + ")")
                    outer = BaseDyadic(x.args[0], y.args[1])
                    outdyad += vect_dot * self.components[x] * \
                               other.components[y] * outer
            return outdyad
        else:
            raise TypeError(str(type(other)) + " not supported for & with dyadics")

    dot = __and__

    @call_highest_priority('__and__')
    def __rand__(self, other):
        if isinstance(other, Symbol):
            outvec = 0
            for x in self.components:
                vect_dot = Symbol("(" + str(x.args[0]) + "." + str(other) + ")")
                outvec += vect_dot * self.components[x] * x.args[1]
            return outvec
        else:
            raise TypeError(str(type(other)) + " not supported for & with dyadic")


class BaseDyadic(Dyadic):
    """ Class to denote a basic dyadic tensor component """
    def __new__(cls, vector1, vector2):
        if not isinstance(vector1, Symbol) or not isinstance(vector2, Symbol):
            raise TypeError("BaseDyadic cannot be composed of non-base "+ \
                            "vectors")
        elif vector1 == 0 or vector2 == 0:
            return S(0)
        obj = super(BaseDyadic, cls).__new__(cls, vector1, vector2)
        obj._base_dyad = obj
        obj._measure_number = 1
        return obj

    @property
    def components(self):
        return {self : 1}

    def __str__(self, printer=None):
        return "(" + str(self.args[0]) + "|" + str(self.args[1]) + ")"

    _sympystr = __str__
    _sympyrepr = _sympystr


class DyadicMul(Mul, Dyadic):
    """ Products of scalars and BaseDyadics """

    def __new__(cls, *args, **options):
        count = 0
        measure_number = 1
        for x in args:
            if x == 0:
                return S(0)
            elif isinstance(x, BaseDyadic) or isinstance(x, DyadicMul):
                count += 1
                dyad = x
                measure_number *= x._measure_number
            elif isinstance(x, DyadicAdd):
                count += 1
                dyad = x
            else:
                measure_number *= x
        if count > 1:
            raise ValueError("Cannot multiple two dyadics")
        elif count == 0:
            raise ValueError("No dyadics supplied")
        if isinstance(dyad, DyadicAdd):
            newargs = [DyadicMul(measure_number, x) for x in dyad.args]
            return DyadicAdd(*newargs)
        obj = super(DyadicMul, cls).__new__(cls, *args, **options)
        obj._base_dyad = dyad._base_dyad
        obj._measure_number = measure_number
        return obj

    def as_coeff_Mul(self, rational=False):
        return (S(1), self)

    __init__ = Mul.__init__

    @property
    def components(self):
        return {self._base_dyad : self._measure_number}


class DyadicAdd(Dyadic, Add):
    """ Class to hold dyadic sums """

    def __new__(cls, *args, **options):
        lookup = {}
        for x in args:
            if x == 0:
                continue
            for y in x.components:
                if y not in lookup:
                    lookup[y] = 0
                lookup[y] += x.components[y]
        newargs = [DyadicMul(x, lookup[x]) for x in lookup]
        obj = super(DyadicAdd, cls).__new__(cls, *newargs, **options)
        if isinstance(obj, Mul):
            return DyadicMul(*obj.args)
        obj._components = lookup
        return obj

    __init__ = Add.__init__

    @property
    def components(self):
        return self._components


def _dyad_div(one, other):
    """ Helper for division involving dyadics """
    if isinstance(one, Dyadic) and isinstance(other, Dyadic):
        raise TypeError("Cannot divide two dyadics")
    elif isinstance(one, Dyadic):
        return DyadicMul(one, Pow(other, S.NegativeOne))
    else:
        raise TypeError("Cannot divide by a dyadic")

