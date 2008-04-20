""" Defines default methods for unary and binary operations.
"""

from basic import Basic, S
from sympify import sympify, SympifyError, _sympifyit
#from add import Add    /cyclic/
#from mul import Mul    /cyclic/
#from pow import Pow    /cyclic/
#from relational import Equality, Unequality, Inequality, StrictInequality /cyclic/
#from sympy.functions.elementary.complexes import abs as abs_   /cyclic/

def _no_unary_operation(op, obj):
    return 'unary operation `%s` not defined for %s' % (op, obj.__class__.__name__)
def _no_binary_operation(op, obj1, obj2):
    return 'binary operation `%s` not defined between %s and %s' \
           % (op, obj1.__class__.__name__, obj2.__class__.__name__)

class ArithMeths(object):

    __slots__ = []

    def __pos__(self):
        return self
    def __neg__(self):
        return Mul(S.NegativeOne, self)
    def __abs__(self):
        return abs_(self)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        return Add(self, other)
    @_sympifyit('other', NotImplemented)
    def __radd__(self, other):
        return Add(other, self)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        return Add(self, -other)
    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return Add(other, -self)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        return Mul(self, other)
    @_sympifyit('other', NotImplemented)
    def __rmul__(self, other):
        return Mul(other, self)

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        return Pow(self, other)
    @_sympifyit('other', NotImplemented)
    def __rpow__(self, other):
        return Pow(other, self)

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        return Mul(self, Pow(other, S.NegativeOne))
    @_sympifyit('other', NotImplemented)
    def __rdiv__(self, other):
        return Mul(other, Pow(self, S.NegativeOne))

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

class NoArithMeths(object):

    __slots__ = []

    def __pos__(self):
        raise TypeError, _no_unary_operation('+', self)
    def __neg__(self):
        raise TypeError, _no_unary_operation('-', self)
    def __add__(self, other):
        raise TypeError, _no_binary_operation('+', self, other)
    def __radd__(self, other):
        raise TypeError, _no_binary_operation('+', other, self)
    def __sub__(self, other):
        raise TypeError, _no_binary_operation('-', self, other)
    def __rsub__(self, other):
        raise TypeError, _no_binary_operation('-', other, self)
    def __mul__(self, other):
        raise TypeError, _no_binary_operation('*', self, other)
    def __rmul__(self, other):
        raise TypeError, _no_binary_operation('*', other, self)
    def __div__(self, other):
        raise TypeError, _no_binary_operation('/', self, other)
    def __rdiv__(self, other):
        raise TypeError, _no_binary_operation('/', other, self)
    def __pow__(self, other):
        raise TypeError, _no_binary_operation('**', self, other)
    def __rpow__(self, other):
        raise TypeError, _no_binary_operation('**', other, self)
    def _eval_power(self, other):
        return None

class RelMeths(object):

    __slots__ = []

    # TODO all comparison methods should return True/False directly
    # see #153

    @_sympifyit('other', False) # sympy >  other
    def __lt__(self, other):
        #return sympify(other) > self
        return StrictInequality(self, other)

    @_sympifyit('other', True)  # sympy >  other
    def __gt__(self, other):
        return StrictInequality(other, self)
        #return sympify(other) < self

    @_sympifyit('other', False) # sympy >  other
    def __le__(self, other):
        return Inequality(self, other)

    @_sympifyit('other', True)  # sympy >  other
    def __ge__(self, other):
        return sympify(other) <= self

class NoRelMeths(object):

    __slots__ = []

    def __lt__(self, other):
        return hash(self) < hash(other)
        raise TypeError, _no_binary_operation('<', self, other)
    def __gt__(self, other):
        return hash(self) > hash(other)
        raise TypeError, _no_binary_operation('>', self, other)
    def __le__(self, other):
        raise TypeError, _no_binary_operation('<=', self, other)
    def __ge__(self, other):
        raise TypeError, _no_binary_operation('>=', self, other)
