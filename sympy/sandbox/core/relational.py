
from utils import memoizer_immutable_args
from basic import Composite

class Relational(Composite, tuple):

    @memoizer_immutable_args('Relational.__new__')
    def __new__(cls, lhs, rhs):
        sympify = cls.sympify
        lhs, rhs = sympify(lhs), sympify(rhs)
        return tuple.__new__(cls, (lhs, rhs))

    @property
    def lhs(self):
        return self[0]

    @property
    def rhs(self):
        return self[1]

class Equality(Relational):

    rel_op = '=='

    @memoizer_immutable_args('Equality.__new__')
    def __new__(cls, lhs, rhs):
        sympify = cls.sympify
        lhs, rhs = sympify(lhs), sympify(rhs)
        if lhs.compare(rhs)==0: return True
        return tuple.__new__(cls, (lhs, rhs))

    @memoizer_immutable_args('Equality.__nonzero__')
    def __nonzero__(self):
        return self.lhs.compare(self.rhs)==0

class Unequality(Relational):

    rel_op = '!='

    @memoizer_immutable_args('Unequality.__new__')
    def __new__(cls, lhs, rhs):
        sympify = cls.sympify
        lhs, rhs = sympify(lhs), sympify(rhs)
        if lhs.compare(rhs)==0: return False
        return tuple.__new__(cls, (lhs, rhs))

    @memoizer_immutable_args('Unequality.__nonzero__')
    def __nonzero__(self):
        return self.lhs.compare(self.rhs)!=0

class StrictInequality(Relational):

    rel_op = '<'

    @memoizer_immutable_args('StrictInequality.__nonzero__')
    def __nonzero__(self):
        return self.lhs.compare(self.rhs)==-1

class Inequality(Relational):

    rel_op = '<='

    @memoizer_immutable_args('Inequality.__nonzero__')
    def __nonzero__(self):
        return self.lhs.compare(self.rhs)<=0
