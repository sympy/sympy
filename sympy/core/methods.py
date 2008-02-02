""" Defines default methods for unary and binary operations.
"""

from basic import Basic, S, SympifyError, sympify
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

    def __pos__(self):
        return self
    def __neg__(self):
        return S.NegativeOne * self
    def __abs__(self):
        return abs_(self)
    def __add__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        if not isinstance(other, Basic):
            return NotImplemented
        return Add(self, other)
    def __radd__(self, other):
        return sympify(other).__add__(self)
    def __sub__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        if not isinstance(other, Basic):
            return NotImplemented
        return self + (-other)
    def __rsub__(self, other):
        return sympify(other).__sub__(self)
    def __mul__(self, other):
        # FIXME this is a dirty hack. matrix should be ordinary SymPy object
        from sympy.matrices import Matrix
        if isinstance(other, Matrix): return NotImplemented
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        if not isinstance(other, Basic):
            return NotImplemented
        return Mul(self, other)
    def __rmul__(self, other):
        return sympify(other).__mul__(self)
    def __pow__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        if not isinstance(other, Basic):
            return NotImplemented
        return Pow(self, other)
    def __rpow__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        if not isinstance(other, Basic):
            return NotImplemented
        return other.__pow__(self)
    def __div__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        if not isinstance(other, Basic):
            return NotImplemented
        return self * (other ** S.NegativeOne)
    def __truediv__(self, other):
        return self.__div__(other)
    def __rdiv__(self, other):
        return sympify(other).__div__(self)
    def __rtruediv__(self, other):
        return self.__rdiv__(other)

class NoArithMeths(object):
    
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
    
    def __eq__(self, other):
        try:
            other = sympify(other)
        except ValueError:
            return False
        return Equality(self, other)
    def __ne__(self, other):
        try:
            other = sympify(other)
        except ValueError:
            return True
        return Unequality(self, other)
    def __lt__(self, other):
        #return sympify(other) > self
        return StrictInequality(self, other)
    def __gt__(self, other):
        return StrictInequality(other, self)
        #return sympify(other) < self
    def __le__(self, other):
        return Inequality(self, other)
    def __ge__(self, other):
        return sympify(other) <= self

class NoRelMeths(object):

    def __eq__(self, other):
        return Equality(self, other)
    def __ne__(self, other):
        return Unequality(self, other)
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
