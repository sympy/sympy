
from basic import Composite

class Relational(Composite, tuple):

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

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)==0

class Unequality(Relational):

    rel_op = '!='

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)!=0

class StrictInequality(Relational):

    rel_op = '<'

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)==-1

class Inequality(Relational):

    rel_op = '<='

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)<=0
