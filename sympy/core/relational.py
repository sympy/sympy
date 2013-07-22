from .basic import S
from .expr import Expr
from .evalf import EvalfMixin
from .symbol import Symbol
from .sympify import _sympify

from sympy.logic.boolalg import Boolean

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
    return Relational(a, b, op)


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
    return Relational(a, b, '==')


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
    return Relational(a, b, '!=')


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
    return Relational(a, b, '<')


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
    return Relational(a, b, '<=')


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
    return Relational(a, b, '>')


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
    return Relational(a, b, '>=')

# Note, see issue 1887.  Ideally, we wouldn't want to subclass both Boolean
# and Expr.


class Relational(Boolean, Expr, EvalfMixin):

    __slots__ = []

    is_Relational = True

    # ValidRelationOperator - Defined below, because the necessary classes
    #   have not yet been defined

    def __new__(cls, lhs, rhs, rop=None, **assumptions):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        if cls is not Relational:
            rop_cls = cls
        else:
            try:
                rop_cls = cls.ValidRelationOperator[ rop ]
            except KeyError:
                msg = "Invalid relational operator symbol: '%r'"
                raise ValueError(msg % repr(rop))

        diff = S.NaN
        if isinstance(lhs, Expr) and isinstance(rhs, Expr):
            diff = lhs - rhs
        if not (diff is S.NaN or diff.has(Symbol)):
            know = diff.equals(0, failing_expression=True)
            if know is True:  # exclude failing expression case
                diff = S.Zero
            elif know is False:
                diff = diff.n()
        if rop_cls is Equality:
            if (lhs == rhs) is True or (diff == S.Zero) is True:
                return True
            elif diff is S.NaN:
                pass
            elif diff.is_Number or diff.is_Float:
                return False
            elif lhs.is_real is not rhs.is_real and \
                lhs.is_real is not None and \
                   rhs.is_real is not None:
                return False
        elif rop_cls is Unequality:
            if (lhs == rhs) is True or (diff == S.Zero) is True:
                return False
            elif diff is S.NaN:
                pass
            elif diff.is_Number or diff.is_Float:
                return True
            elif lhs.is_real is not rhs.is_real and \
                lhs.is_real is not None and \
                   rhs.is_real is not None:
                return True
        elif diff.is_Number and diff.is_real:
            return rop_cls._eval_relation(diff, S.Zero)

        obj = Expr.__new__(rop_cls, lhs, rhs, **assumptions)
        return obj

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

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

    def _eval_simplify(self, ratio, measure):
        return self.__class__(self.lhs.simplify(ratio=ratio),
                              self.rhs.simplify(ratio=ratio))


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
        return self.lhs.compare(self.rhs) == 0

    __bool__ = __nonzero__


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
        return self.lhs.compare(self.rhs) != 0

    __bool__ = __nonzero__


class _Greater(Relational):
    """Not intended for general use

    _Greater is only used so that GreaterThan and StrictGreaterThan may subclass
    it for the .gts and .lts properties.
    """

    __slots__ = ()

    @property
    def gts(self):
        return self._args[0]

    @property
    def lts(self):
        return self._args[1]


class _Less(Relational):
    """Not intended for general use.

    _Less is only used so that LessThan and StrictLessThan may subclass it for
    the .gts and .lts properties.
    """

    __slots__ = ()

    @property
    def gts(self):
        return self._args[1]

    @property
    def lts(self):
        return self._args[0]


class GreaterThan(_Greater):
    """Class representations of inequalities.

    Extended Summary
    ================

    The ``*Than`` classes represent inequal relationships, where the left-hand
    side is generally bigger or smaller than the right-hand side.  For example,
    the GreaterThan class represents an inequal relationship where the
    left-hand side is at least as big as the right side, if not bigger.  In
    mathematical notation:

    lhs >= rhs

    In total, there are four ``*Than`` classes, to represent the four
    inequalities:

    +-----------------+--------+
    |Class Name       | Symbol |
    +=================+========+
    |GreaterThan      | (>=)   |
    +-----------------+--------+
    |LessThan         | (<=)   |
    +-----------------+--------+
    |StrictGreaterThan| (>)    |
    +-----------------+--------+
    |StrictLessThan   | (<)    |
    +-----------------+--------+

    All classes take two arguments, lhs and rhs.

    +----------------------------+-----------------+
    |Signature Example           | Math equivalent |
    +============================+=================+
    |GreaterThan(lhs, rhs)       |   lhs >= rhs    |
    +----------------------------+-----------------+
    |LessThan(lhs, rhs)          |   lhs <= rhs    |
    +----------------------------+-----------------+
    |StrictGreaterThan(lhs, rhs) |   lhs >  rhs    |
    +----------------------------+-----------------+
    |StrictLessThan(lhs, rhs)    |   lhs <  rhs    |
    +----------------------------+-----------------+

    In addition to the normal .lhs and .rhs of Relations, ``*Than`` inequality
    objects also have the .lts and .gts properties, which represent the "less
    than side" and "greater than side" of the operator.  Use of .lts and .gts
    in an algorithm rather than .lhs and .rhs as an assumption of inequality
    direction will make more explicit the intent of a certain section of code,
    and will make it similarly more robust to client code changes:

    >>> from sympy import GreaterThan, StrictGreaterThan
    >>> from sympy import LessThan,    StrictLessThan
    >>> from sympy import And, Ge, Gt, Le, Lt, Rel, S
    >>> from sympy.abc import x, y, z
    >>> from sympy.core.relational import Relational

    >>> e = GreaterThan(x, 1)
    >>> e
    x >= 1
    >>> '%s >= %s is the same as %s <= %s' % (e.gts, e.lts, e.lts, e.gts)
    'x >= 1 is the same as 1 <= x'

    Examples
    ========

    One generally does not instantiate these classes directly, but uses various
    convenience methods:

    >>> e1 = Ge( x, 2 )      # Ge is a convenience wrapper
    >>> print(e1)
    x >= 2

    >>> rels = Ge( x, 2 ), Gt( x, 2 ), Le( x, 2 ), Lt( x, 2 )
    >>> print('%s\\n%s\\n%s\\n%s' % rels)
    x >= 2
    x > 2
    x <= 2
    x < 2

    Another option is to use the Python inequality operators (>=, >, <=, <)
    directly.  Their main advantage over the Ge, Gt, Le, and Lt counterparts, is
    that one can write a more "mathematical looking" statement rather than
    littering the math with oddball function calls.  However there are certain
    (minor) caveats of which to be aware (search for 'gotcha', below).

    >>> e2 = x >= 2
    >>> print(e2)
    x >= 2
    >>> print("e1: %s,    e2: %s" % (e1, e2))
    e1: x >= 2,    e2: x >= 2
    >>> e1 == e2
    True

    However, it is also perfectly valid to instantiate a ``*Than`` class less
    succinctly and less conveniently:

    >>> rels = Rel(x, 1, '>='), Relational(x, 1, '>='), GreaterThan(x, 1)
    >>> print('%s\\n%s\\n%s' % rels)
    x >= 1
    x >= 1
    x >= 1

    >>> rels = Rel(x, 1, '>'), Relational(x, 1, '>'), StrictGreaterThan(x, 1)
    >>> print('%s\\n%s\\n%s' % rels)
    x > 1
    x > 1
    x > 1

    >>> rels = Rel(x, 1, '<='), Relational(x, 1, '<='), LessThan(x, 1)
    >>> print("%s\\n%s\\n%s" % rels)
    x <= 1
    x <= 1
    x <= 1

    >>> rels = Rel(x, 1, '<'), Relational(x, 1, '<'), StrictLessThan(x, 1)
    >>> print('%s\\n%s\\n%s' % rels)
    x < 1
    x < 1
    x < 1

    Notes
    =====

    There are a couple of "gotchas" when using Python's operators.

    The first enters the mix when comparing against a literal number as the lhs
    argument.  Due to the order that Python decides to parse a statement, it may
    not immediately find two objects comparable.  For example, to evaluate the
    statement (1 < x), Python will first recognize the number 1 as a native
    number, and then that x is *not* a native number.  At this point, because a
    native Python number does not know how to compare itself with a SymPy object
    Python will try the reflective operation, (x > 1).  Unfortunately, there is
    no way available to SymPy to recognize this has happened, so the statement
    (1 < x) will turn silently into (x > 1).

    >>> e1 = x >  1
    >>> e2 = x >= 1
    >>> e3 = x <  1
    >>> e4 = x <= 1
    >>> e5 = 1 >  x
    >>> e6 = 1 >= x
    >>> e7 = 1 <  x
    >>> e8 = 1 <= x
    >>> print("%s     %s\\n"*4 % (e1, e2, e3, e4, e5, e6, e7, e8))
    x > 1     x >= 1
    x < 1     x <= 1
    x < 1     x <= 1
    x > 1     x >= 1

    If the order of the statement is important (for visual output to the
    console, perhaps), one can work around this annoyance in a couple ways: (1)
    "sympify" the literal before comparison, (2) use one of the wrappers, or (3)
    use the less succinct methods described above:

    >>> e1 = S(1) >  x
    >>> e2 = S(1) >= x
    >>> e3 = S(1) <  x
    >>> e4 = S(1) <= x
    >>> e5 = Gt(1, x)
    >>> e6 = Ge(1, x)
    >>> e7 = Lt(1, x)
    >>> e8 = Le(1, x)
    >>> print("%s     %s\\n"*4 % (e1, e2, e3, e4, e5, e6, e7, e8))
    1 > x     1 >= x
    1 < x     1 <= x
    1 > x     1 >= x
    1 < x     1 <= x

    The other gotcha is with chained inequalities.  Occasionally, one may be
    tempted to write statements like:

    >>> e = x < y < z  # silent error!  Where did ``x`` go?
    >>> e #doctest: +SKIP
    y < z

    Due to an implementation detail or decision of Python [1]_, there is no way
    for SymPy to reliably create that as a chained inequality.  To create a
    chained inequality, the only method currently available is to make use of
    And:

    >>> e = And(x < y, y < z)
    >>> type( e )
    And
    >>> e
    And(x < y, y < z)

    Note that this is different than chaining an equality directly via use of
    parenthesis (this is currently an open bug in SymPy [2]_):

    >>> e = (x < y) < z
    >>> type( e )
    <class 'sympy.core.relational.StrictLessThan'>
    >>> e
    (x < y) < z

    Any code that explicitly relies on this latter functionality will not be
    robust as this behaviour is completely wrong and will be corrected at some
    point.  For the time being (circa Jan 2012), use And to create chained
    inequalities.

    .. [1] This implementation detail is that Python provides no reliable
       method to determine that a chained inequality is being built.  Chained
       comparison operators are evaluated pairwise, using "and" logic (see
       http://docs.python.org/reference/expressions.html#notin).  This is done
       in an efficient way, so that each object being compared is only
       evaluated once and the comparison can short-circuit.  For example, ``1
       > 2 > 3`` is evaluated by Python as ``(1 > 2) and (2 > 3)``.  The
       ``and`` operator coerces each side into a bool, returning the object
       itself when it short-circuits.  Currently, the bool of the --Than
       operators will give True or False arbitrarily.  Thus, if we were to
       compute ``x > y > z``, with ``x``, ``y``, and ``z`` being Symbols,
       Python converts the statement (roughly) into these steps:

        (1) x > y > z
        (2) (x > y) and (y > z)
        (3) (GreaterThanObject) and (y > z)
        (4) (GreaterThanObject.__nonzero__()) and (y > z)
        (5) (True) and (y > z)
        (6) (y > z)
        (7) LessThanObject

       Because of the "and" added at step 2, the statement gets turned into a
       weak ternary statement.  If the first object evalutes __nonzero__ as
       True, then the second object, (y > z) is returned.  If the first object
       evaluates __nonzero__ as False (step 5), then (x > y) is returned.

           In Python, there is no way to override the ``and`` operator, or to
           control how it short circuits, so it is impossible to make something
           like ``x > y > z`` work.  There is an open PEP to change this,
           :pep:`335`, but until that is implemented, this cannot be made to work.

    .. [2] For more information, see these two bug reports:

       "Separate boolean and symbolic relationals"
       `Issue 1887 <http://code.google.com/p/sympy/issues/detail?id=1887>`_

       "It right 0 < x < 1 ?"
       `Issue 2960 <http://code.google.com/p/sympy/issues/detail?id=2960>`_

    """

    rel_op = '>='

    __slots__ = ()

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs >= rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) >= 0

    __bool__ = __nonzero__


class LessThan(_Less):
    __doc__ = GreaterThan.__doc__
    __slots__ = ()

    rel_op = '<='

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs <= rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) <= 0

    __bool__ = __nonzero__


class StrictGreaterThan(_Greater):
    __doc__ = GreaterThan.__doc__
    __slots__ = ()

    rel_op = '>'

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs > rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) > 0

    __bool__ = __nonzero__


class StrictLessThan(_Less):
    __doc__ = GreaterThan.__doc__
    __slots__ = ()

    rel_op = '<'

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return lhs < rhs

    def __nonzero__(self):
        return self.lhs.compare( self.rhs ) < 0

    __bool__ = __nonzero__

# A class-specific (not object-specific) data item used for a minor speedup.  It
# is defined here, rather than directly in the class, because the classes that
# it references have not been defined until now (e.g. StrictLessThan).
Relational.ValidRelationOperator = {
    None: Equality,
    '==': Equality,
    'eq': Equality,
    '!=': Unequality,
    '<>': Unequality,
    'ne': Unequality,
    '>=': GreaterThan,
    'ge': GreaterThan,
    '<=': LessThan,
    'le': LessThan,
    '>': StrictGreaterThan,
    'gt': StrictGreaterThan,
    '<': StrictLessThan,
    'lt': StrictLessThan,
}
