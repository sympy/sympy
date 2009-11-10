
from basic import Basic
from sympify import _sympify

def Rel(a, b, op):
    """
    A handy wrapper around the Relational class.
    Rel(a,b, op)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Rel(y, x+x**2, '==')
    y == x + x**2

    """
    return Relational(a,b,op)

def Eq(a, b=0):
    """
    A handy wrapper around the Relational class.
    Eq(a,b)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Eq(y, x+x**2)
    y == x + x**2

    """
    return Relational(a,b,'==')

def Ne(a, b):
    """
    A handy wrapper around the Relational class.
    Ne(a,b)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Ne(y, x+x**2)
    y != x + x**2

    """
    return Relational(a,b,'!=')

def Lt(a, b):
    """
    A handy wrapper around the Relational class.
    Lt(a,b)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Lt(y, x+x**2)
    y < x + x**2

    """
    return Relational(a,b,'<')

def Le(a, b):
    """
    A handy wrapper around the Relational class.
    Le(a,b)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Le(y, x+x**2)
    y <= x + x**2

    """
    return Relational(a,b,'<=')

def Gt(a, b):
    """
    A handy wrapper around the Relational class.
    Gt(a,b)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Gt(y, x+x**2)
    x + x**2 < y

    """
    return Relational(a,b,'>')

def Ge(a, b):
    """
    A handy wrapper around the Relational class.
    Ge(a,b)

    Example:
    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> Ge(y, x+x**2)
    x + x**2 <= y

    """
    return Relational(a,b,'>=')

class Relational(Basic):

    __slots__ = []

    @staticmethod
    def get_relational_class(rop):
        if rop is None or rop in ['==','eq']: return Equality, False
        if rop in ['!=','<>','ne']: return Unequality, False
        if rop in ['<','lt']: return StrictInequality, False
        if rop in ['>','gt']: return StrictInequality, True
        if rop in ['<=','le']: return Inequality, False
        if rop in ['>=','ge']: return Inequality, True
        raise ValueError("Invalid relational operator symbol: %r" % (rop))

    def __new__(cls, lhs, rhs, rop=None, **assumptions):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
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

    def _eval_subs(self, old, new):
        return self.__class__(self.lhs._eval_subs(old, new), self.rhs._eval_subs(old, new))

class Equality(Relational):

    rel_op = '=='

    __slots__ = []

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)==0

class Unequality(Relational):

    rel_op = '!='

    __slots__ = []

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)!=0

class StrictInequality(Relational):

    rel_op = '<'

    __slots__ = []

    def __nonzero__(self):
        if self.lhs.is_comparable and self.rhs.is_comparable:
            if self.lhs.is_Number and self.rhs.is_Number:
                return self.lhs < self.rhs
            return self.lhs.evalf()<self.rhs.evalf()
        return self.lhs.compare(self.rhs)==-1

class Inequality(Relational):

    rel_op = '<='

    __slots__ = []

    def __nonzero__(self):
        if self.lhs.is_comparable and self.rhs.is_comparable:
            if self.lhs.is_Number and self.rhs.is_Number:
                return self.lhs <= self.rhs
            return self.lhs.evalf()<=self.rhs.evalf()
        return self.lhs.compare(self.rhs)<=0
