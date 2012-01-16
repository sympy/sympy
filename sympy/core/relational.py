from expr import Expr
from evalf import EvalfMixin
from sympify import _sympify

__all__ = (
 'Rel', 'Eq', 'Ne', 'Lt', 'Le', 'Gt', 'Ge',
 'Relational', 'Equality', 'Unequality', 'StrictLessThan', 'LessThan',
 'StrictGreaterThan', 'GreaterThan',
)


def Rel(a, b, op):
    """
    A handy wrapper around the Relational class.
    Rel(a,b, op)

    Examples
    ========

    >>> from sympy import Rel
    >>> from sympy.abc import x, y
    >>> Rel(y, x+x**2, '==')
    y == x**2 + x

    """
    return Relational(a,b,op)

def Eq(a, b=0):
    """
    A handy wrapper around the Relational class.
    Eq(a,b)

    Examples
    ========

    >>> from sympy import Eq
    >>> from sympy.abc import x, y
    >>> Eq(y, x+x**2)
    y == x**2 + x

    """
    return Relational(a,b,'==')

def Ne(a, b):
    """
    A handy wrapper around the Relational class.
    Ne(a,b)

    Examples
    ========

    >>> from sympy import Ne
    >>> from sympy.abc import x, y
    >>> Ne(y, x+x**2)
    y != x**2 + x

    """
    return Relational(a,b,'!=')

def Lt(a, b):
    """
    A handy wrapper around the Relational class.
    Lt(a,b)

    Examples
    ========

    >>> from sympy import Lt
    >>> from sympy.abc import x, y
    >>> Lt(y, x+x**2)
    y < x**2 + x

    """
    return Relational(a,b,'<')

def Le(a, b):
    """
    A handy wrapper around the Relational class.
    Le(a,b)

    Examples
    ========

    >>> from sympy import Le
    >>> from sympy.abc import x, y
    >>> Le(y, x+x**2)
    y <= x**2 + x

    """
    return Relational(a,b,'<=')

def Gt(a, b):
    """
    A handy wrapper around the Relational class.
    Gt(a,b)

    Examples
    ========

    >>> from sympy import Gt
    >>> from sympy.abc import x, y
    >>> Gt(y, x + x**2)
    y > x**2 + x

    """
    return Relational(a,b,'>')

def Ge(a, b):
    """
    A handy wrapper around the Relational class.
    Ge(a,b)

    Examples
    ========

    >>> from sympy import Ge
    >>> from sympy.abc import x, y
    >>> Ge(y, x + x**2)
    y >= x**2 + x

    """
    return Relational(a,b,'>=')

class Relational(Expr, EvalfMixin):

    __slots__ = []

    is_Relational = True

    # ValidRelationOperator - Defined below, because the necessary classes
    #   have not yet been defined

    @staticmethod
    def get_relational_class(rop):
        try:
            return Relational.ValidRelationOperator[ rop ]
        except KeyError:
            raise ValueError("Invalid relational operator symbol: %r" % (rop))

    def __new__(cls, lhs, rhs, rop=None, **assumptions):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        if cls is not Relational:
            rop_cls = cls
        else:
            rop_cls = Relational.get_relational_class(rop)
        if lhs.is_number and rhs.is_number and lhs.is_real and rhs.is_real:
            # Just becase something is a number, doesn't mean you can evalf it.
            Nlhs = lhs.evalf()
            if Nlhs.is_Number:
                # S.Zero.evalf() returns S.Zero, so test Number instead of Float
                Nrhs = rhs.evalf()
                if Nrhs.is_Number:
                    return rop_cls._eval_relation(Nlhs, Nrhs)

        obj = Expr.__new__(rop_cls, lhs, rhs, **assumptions)
        return obj

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    def _eval_subs(self, old, new):
        if self == old:
            return new
        return self.__class__(self.lhs._eval_subs(old, new), self.rhs._eval_subs(old, new))

    def _eval_evalf(self, prec):
        return self.func(*[s._evalf(prec) for s in self.args])

    def doit(self, **hints):
        lhs = self.lhs
        rhs = self.rhs
        if hints.get('deep', True):
            lhs = lhs.doit(**hints)
            rhs = rhs.doit(**hints)
        return self._eval_relation_doit(lhs, rhs)

    @classmethod
    def _eval_relation_doit(cls, lhs, rhs):
        return cls._eval_relation(lhs, rhs)

class Equality(Relational):

    rel_op = '=='

    __slots__ = []

    is_Equality = True

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs == rhs

    @classmethod
    def _eval_relation_doit(cls, lhs, rhs):
        return Eq(lhs, rhs)

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)==0

class Unequality(Relational):

    rel_op = '!='

    __slots__ = []

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs != rhs

    @classmethod
    def _eval_relation_doit(cls, lhs, rhs):
        return Ne(lhs, rhs)

    def __nonzero__(self):
        return self.lhs.compare(self.rhs)!=0

class _Greater(Relational):
    @property
    def gts(self):
        return self._args[0]

    @property
    def lts(self):
        return self._args[1]

class _Less(Relational):
    @property
    def gts(self):
        return self._args[1]

    @property
    def lts(self):
        return self._args[0]

    def __eq__ ( self, other ):
        if isinstance(other, _Less):
            ot = other._hashable_content()
        elif isinstance(other, _Greater):
            ot = tuple(reversed( other._hashable_content() ))
        else:
            return False

        st = self._hashable_content()

        return st == ot

class GreaterThan(_Greater):

    rel_op = '>='

    __slots__ = []

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs >= rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) >= 0

    def __eq__ ( self, other ):
        if isinstance(other, GreaterThan):
            ot = other._hashable_content()
        elif isinstance(other, LessThan):
            ot = tuple(reversed( other._hashable_content() ))
        else:
            return False

        st = self._hashable_content()

        return st == ot

class LessThan(_Less):

    rel_op = '<='

    __slots__ = []

    @classmethod
    def _eval_relation(cls, lhs, rhs):
         return lhs <= rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) <= 0

    def __eq__ ( self, other ):
        if isinstance(other, LessThan):
            ot = other._hashable_content()
        elif isinstance(other, GreaterThan):
            ot = tuple(reversed( other._hashable_content() ))
        else:
            return False

        st = self._hashable_content()

        return st == ot

class StrictGreaterThan(_Greater):

    rel_op = '>'

    __slots__ = []

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs > rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) == 1

    def __eq__ ( self, other ):
        if isinstance(other, StrictGreaterThan):
            ot = other._hashable_content()
        elif isinstance(other, StrictLessThan):
            ot = tuple(reversed( other._hashable_content() ))
        else:
            return False

        st = self._hashable_content()

        return st == ot

class StrictLessThan(_Less):

    rel_op = '<'

    __slots__ = []

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs < rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) == -1

    def __eq__ ( self, other ):
        if isinstance(other, StrictLessThan):
            ot = other._hashable_content()
        elif isinstance(other, StrictGreaterThan):
            ot = tuple(reversed( other._hashable_content() ))
        else:
            return False

        st = self._hashable_content()

        return st == ot

# A class (not object) specific data item used for a minor speedup.  It is
# defined here, rather than directly in the class, because the classes that it
# references have not been defined until now (e.g. StrictLessThan).
Relational.ValidRelationOperator = {
  None : Equality,
  '==' : Equality,
  'eq' : Equality,
  '!=' : Unequality,
  '<>' : Unequality,
  'ne' : Unequality,
  '>=' : GreaterThan,
  'ge' : GreaterThan,
  '<=' : LessThan,
  'le' : LessThan,
  '>'  : StrictGreaterThan,
  'gt' : StrictGreaterThan,
  '<'  : StrictLessThan,
  'lt' : StrictLessThan,
}
