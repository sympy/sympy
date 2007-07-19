
from basic import Basic
from methods import NoRelMeths

class Relational(Basic, NoRelMeths):

    @staticmethod
    def get_relational_class(rop):
        if rop is None or rop in ['==','eq']: return Equality, False
        if rop in ['!=','<>','ne']: return UnEquality, False
        if rop in ['<','lt']: return StrictInequality, False
        if rop in ['>','gt']: return StrictInequality, True
        if rop in ['<=','le']: return Inequality, False
        if rop in ['>=','ge']: return Inequality, True
        raise ValueError("Invalid relational operator symbol: %r" % (rop))

    def __new__(cls, lhs, rhs, rop=None, **assumptions):
        lhs = Basic.sympify(lhs)
        rhs = Basic.sympify(rhs)
        if cls is not Relational:
            rop_cls = cls
        else:
            rop_cls, swap = Relational.get_relational_class(rop)
            if swap: lhs, rhs = rhs, lhs
        obj = Basic.__new__(rop_cls, lhs, rhs, **assumptions)
        return obj

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    @property
    def precedence(self):
        return 20

    def tostr(self, level=0):
        precedence = self.precedence
        r = '%s %s %s' % (self.lhs.tostr(precedence), self.rel_op, self.rhs.tostr(precedence))
        if precedence <= level:
            return '(%s)' % r
        return r

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
        if self.lhs.is_comparable and self.rhs.is_comparable:
            if isinstance(self.lhs, Basic.Number) and isinstance(self.rhs, Basic.Number):
                return self.lhs < self.rhs
            return self.lhs.evalf()<self.rhs.evalf()
        return self.lhs.compare(self.rhs)==-1

class Inequality(Relational):

    rel_op = '<='

    def __nonzero__(self):
        if self.lhs.is_comparable and self.rhs.is_comparable:
            if isinstance(self.lhs, Basic.Number) and isinstance(self.rhs, Basic.Number):
                return self.lhs <= self.rhs
            return self.lhs.evalf()<=self.rhs.evalf()
        return self.lhs.compare(self.rhs)<=0
