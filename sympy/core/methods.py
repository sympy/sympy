""" Defines default methods for unary and binary operations.
"""

from basic import Basic, cache_it

def _no_unary_operation(op, obj):
    return 'unary operation `%s` not defined for %s' % (op, obj.__class__.__name__)
def _no_binary_operation(op, obj1, obj2):
    return 'binary operation `%s` not defined between %s and %s' \
           % (op, obj1.__class__.__name__, obj2.__class__.__name__)

class ArithMeths(object):

    def __pos__(self):
        return self
    def __neg__(self):
        return Basic.Integer(-1) * self
    def __abs__(self):
        return Basic.Abs()(self)
    def __add__(self, other):
        return Basic.Add(self, other)
    def __radd__(self, other):
        return Basic.sympify(other) + self
    def __sub__(self, other):
        return self + (-Basic.sympify(other))
    def __rsub__(self, other):
        return Basic.sympify(other) - self
    def __mul__(self, other):
        return Basic.Mul(self, other)
    def __rmul__(self, other):
        return Basic.sympify(other) * self
    def __pow__(self, other):
        return Basic.Pow(self, other)
    def __rpow__(self, other):
        return Basic.sympify(other) ** self
    def __div__(self, other):
        return self * (Basic.sympify(other) ** Basic.Integer(-1))
    def __rdiv__(self, other):
        return Basic.sympify(other) / self
    def _eval_power(self, other):
        """ Evaluate Pow(self, other), return new object or None if
        no evaluation can be carried out. This method can be called
        only from Pow.__new__.
        """
        return None

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
        return Basic.Equality(self, other)
    def __ne__(self, other):
        return Basic.Unequality(self, other)
    def __lt__(self, other):
        return Basic.StrictInequality(self, other)
    def __gt__(self, other):
        return Basic.sympify(other) < self
    def __le__(self, other):
        return Basic.Inequality(self, other)
    def __ge__(self, other):
        return Basic.sympify(other) <= self

class NoRelMeths(object):

    def __eq__(self, other):
        return Basic.Equality(self, other)
    def __ne__(self, other):
        return Basic.Unequality(self, other)
    def __lt__(self, other):
        raise TypeError, _no_binary_operation('<', self, other)
    def __gt__(self, other):
        raise TypeError, _no_binary_operation('>', self, other)
    def __le__(self, other):
        raise TypeError, _no_binary_operation('<=', self, other)
    def __ge__(self, other):
        raise TypeError, _no_binary_operation('>=', self, other)
