from __future__ import print_function, division

from .basic import S
from .compatibility import ordered
from .expr import Expr
from .evalf import EvalfMixin
from .function import _coeff_isneg
from .sympify import _sympify
from .evaluate import global_evaluate

from sympy.logic.boolalg import Boolean, BooleanAtom

__all__ = (
    'Rel', 'Eq', 'Ne', 'Lt', 'Le', 'Gt', 'Ge',
    'Relational', 'Equality', 'Unequality', 'StrictLessThan', 'LessThan',
    'StrictGreaterThan', 'GreaterThan',
)


# Note, see issue 4986.  Ideally, we wouldn't want to subclass both Boolean
# and Expr.

class Relational(Boolean, Expr, EvalfMixin):
    """Base class for all relation types.

    Subclasses of Relational should generally be instantiated directly, but
    Relational can be instantiated with a valid `rop` value to dispatch to
    the appropriate subclass.

    Parameters
    ==========
    rop : str or None
        Indicates what subclass to instantiate.  Valid values can be found
        in the keys of Relational.ValidRelationalOperator.

    Examples
    ========

    >>> from sympy import Rel
    >>> from sympy.abc import x, y
    >>> Rel(y, x+x**2, '==')
    Eq(y, x**2 + x)

    """
    __slots__ = []

    is_Relational = True

    # ValidRelationOperator - Defined below, because the necessary classes
    #   have not yet been defined

    def __new__(cls, lhs, rhs, rop=None, **assumptions):
        # If called by a subclass, do nothing special and pass on to Expr.
        if cls is not Relational:
            return Expr.__new__(cls, lhs, rhs, **assumptions)
        # If called directly with an operator, look up the subclass
        # corresponding to that operator and delegate to it
        try:
            cls = cls.ValidRelationOperator[rop]
            return cls(lhs, rhs, **assumptions)
        except KeyError:
            raise ValueError("Invalid relational operator symbol: %r" % rop)

    @property
    def lhs(self):
        """The left-hand side of the relation."""
        return self._args[0]

    @property
    def rhs(self):
        """The right-hand side of the relation."""
        return self._args[1]

    @property
    def reversed(self):
        """Return the relationship with sides (and sign) reversed.

        Examples
        ========

        >>> from sympy import Eq
        >>> from sympy.abc import x
        >>> Eq(x, 1)
        Eq(x, 1)
        >>> _.reversed
        Eq(1, x)
        >>> x < 1
        x < 1
        >>> _.reversed
        1 > x
        """
        ops = {Gt: Lt, Ge: Le, Lt: Gt, Le: Ge}
        a, b = self.args
        return ops.get(self.func, self.func)(b, a, evaluate=False)

    def _eval_evalf(self, prec):
        return self.func(*[s._evalf(prec) for s in self.args])

    @property
    def canonical(self):
        """Return a canonical form of the relational.

        The rules for the canonical form, in order of decreasing priority are:
            1) Number on right if left is not a Number;
            2) Symbol on the left;
            3) Gt/Ge changed to Lt/Le;
            4) Lt/Le are unchanged;
            5) Eq and Ne get ordered args.
        """
        r = self
        if r.func in (Ge, Gt):
            r = r.reversed
        elif r.func in (Lt, Le):
            pass
        elif r.func in (Eq, Ne):
            r = r.func(*ordered(r.args), evaluate=False)
        else:
            raise NotImplemented
        if r.lhs.is_Number and not r.rhs.is_Number:
            r = r.reversed
        elif r.rhs.is_Symbol and not r.lhs.is_Symbol:
            r = r.reversed
        if _coeff_isneg(r.lhs):
            r = r.reversed.func(-r.lhs, -r.rhs, evaluate=False)
        return r

    def equals(self, other, failing_expression=False):
        """Return True if the sides of the relationship are mathematically
        identical and the type of relationship is the same.
        If failing_expression is True, return the expression whose truth value
        was unknown."""
        if isinstance(other, Relational):
            if self == other or self.reversed == other:
                return True
            a, b = self, other
            if a.func in (Eq, Ne) or b.func in (Eq, Ne):
                if a.func != b.func:
                    return False
                l, r = [i.equals(j, failing_expression=failing_expression)
                    for i, j in zip(a.args, b.args)]
                if l is True:
                    return r
                if r is True:
                    return l
                lr, rl = [i.equals(j, failing_expression=failing_expression)
                    for i, j in zip(a.args, b.reversed.args)]
                if lr is True:
                    return rl
                if rl is True:
                    return lr
                e = (l, r, lr, rl)
                if all(i is False for i in e):
                    return False
                for i in e:
                    if i not in (True, False):
                        return i
            else:
                if b.func != a.func:
                    b = b.reversed
                if a.func != b.func:
                    return False
                l = a.lhs.equals(b.lhs, failing_expression=failing_expression)
                if l is False:
                    return False
                r = a.rhs.equals(b.rhs, failing_expression=failing_expression)
                if r is False:
                    return False
                if l is True:
                    return r
                return l

    def _eval_simplify(self, ratio, measure):
        r = self
        r = r.func(*[i.simplify(ratio=ratio, measure=measure)
            for i in r.args])
        if r.is_Relational:
            dif = r.lhs - r.rhs
            # replace dif with a valid Number that will
            # allow a definitive comparison with 0
            v = None
            if dif.is_comparable:
                v = dif.n(2)
            elif dif.equals(0):  # XXX this is expensive
                v = S.Zero
            if v is not None:
                r = r.func._eval_relation(v, S.Zero)

        r = r.canonical
        if measure(r) < ratio*measure(self):
            return r
        else:
            return self

    def __nonzero__(self):
        raise TypeError("cannot determine truth value of Relational")

    __bool__ = __nonzero__

    def as_set(self):
        """
        Rewrites univariate inequality in terms of real sets

        Examples
        ========

        >>> from sympy import Symbol, Eq
        >>> x = Symbol('x', real=True)
        >>> (x>0).as_set()
        (0, oo)
        >>> Eq(x, 0).as_set()
        {0}

        """
        from sympy.solvers.inequalities import solve_univariate_inequality
        syms = self.free_symbols

        if len(syms) == 1:
            sym = syms.pop()
        else:
            raise NotImplementedError("Sorry, Relational.as_set procedure"
                                      " is not yet implemented for"
                                      " multivariate expressions")

        return solve_univariate_inequality(self, sym, relational=False)


Rel = Relational


class Equality(Relational):
    """An equal relation between two objects.

    Represents that two objects are equal.  If they can be easily shown
    to be definitively equal (or unequal), this will reduce to True (or
    False).  Otherwise, the relation is maintained as an unevaluated
    Equality object.  Use the ``simplify`` function on this object for
    more nontrivial evaluation of the equality relation.

    As usual, the keyword argument ``evaluate=False`` can be used to
    prevent any evaluation.

    Examples
    ========

    >>> from sympy import Eq, simplify, exp, cos
    >>> from sympy.abc import x, y
    >>> Eq(y, x + x**2)
    Eq(y, x**2 + x)
    >>> Eq(2, 5)
    False
    >>> Eq(2, 5, evaluate=False)
    Eq(2, 5)
    >>> _.doit()
    False
    >>> Eq(exp(x), exp(x).rewrite(cos))
    Eq(exp(x), sinh(x) + cosh(x))
    >>> simplify(_)
    True

    See Also
    ========

    sympy.logic.boolalg.Equivalent : for representing equality between two
        boolean expressions

    Notes
    =====

    This class is not the same as the == operator.  The == operator tests
    for exact structural equality between two expressions; this class
    compares expressions mathematically.

    If either object defines an `_eval_Eq` method, it can be used in place of
    the default algorithm.  If `lhs._eval_Eq(rhs)` or `rhs._eval_Eq(lhs)`
    returns anything other than None, that return value will be substituted for
    the Equality.  If None is returned by `_eval_Eq`, an Equality object will
    be created as usual.

    """
    rel_op = '=='

    __slots__ = []

    is_Equality = True

    def __new__(cls, lhs, rhs=0, **options):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)

        evaluate = options.pop('evaluate', global_evaluate[0])

        if evaluate:
            # If one expression has an _eval_Eq, return its results.
            if hasattr(lhs, '_eval_Eq'):
                r = lhs._eval_Eq(rhs)
                if r is not None:
                    return r
            if hasattr(rhs, '_eval_Eq'):
                r = rhs._eval_Eq(lhs)
                if r is not None:
                    return r
            # If expressions have the same structure, they must be equal.
            if lhs == rhs:
                return S.true
            elif all(isinstance(i, BooleanAtom) for i in (rhs, lhs)):
                return S.false

            # check finiteness
            fin = L, R = [i.is_finite for i in (lhs, rhs)]
            if None not in fin:
                if L != R:
                    return S.false
                if L is False:
                    return S.true

            # Detect incompatibility such as lhs real and rhs not real.
            if lhs.is_complex and rhs.is_complex:
                r = (lhs - rhs).is_zero
                if r is not None:
                    return _sympify(r)

            # check infiniteness if the difference evaluates
            try:
                n, d = (lhs - rhs).as_numer_denom()
                fin = [i.is_finite for i in (n, d)]
                if None not in fin:
                    if fin[0] != fin[1]:
                        return _sympify(fin[0])
                    elif fin[0]:
                        r = n.is_zero
                        if r is True:
                            if d.is_zero is False:
                                return S.true
                        elif r is False:
                            return S.false
            except (TypeError, AttributeError):
                pass

        return Relational.__new__(cls, lhs, rhs, **options)

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return _sympify(lhs == rhs)

Eq = Equality


class Unequality(Relational):
    """An unequal relation between two objects.

    Represents that two objects are not equal.  If they can be shown to be
    definitively equal, this will reduce to False; if definitively unequal,
    this will reduce to True.  Otherwise, the relation is maintained as an
    Unequality object.

    Examples
    ========

    >>> from sympy import Ne
    >>> from sympy.abc import x, y
    >>> Ne(y, x+x**2)
    Ne(y, x**2 + x)

    See Also
    ========
    Equality

    Notes
    =====
    This class is not the same as the != operator.  The != operator tests
    for exact structural equality between two expressions; this class
    compares expressions mathematically.

    This class is effectively the inverse of Equality.  As such, it uses the
    same algorithms, including any available `_eval_Eq` methods.

    """
    rel_op = '!='

    __slots__ = []

    def __new__(cls, lhs, rhs, **options):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)

        evaluate = options.pop('evaluate', global_evaluate[0])

        if evaluate:
            is_equal = Equality(lhs, rhs)
            if isinstance(is_equal, BooleanAtom):
                return ~is_equal

        return Relational.__new__(cls, lhs, rhs, **options)

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return _sympify(lhs != rhs)

Ne = Unequality


class _Inequality(Relational):
    """Internal base class for all *Than types.

    Each subclass must implement _eval_relation to provide the method for
    comparing two real numbers.

    """
    __slots__ = []

    def __new__(cls, lhs, rhs, **options):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)

        evaluate = options.pop('evaluate', global_evaluate[0])

        if evaluate:
            # First we invoke the appropriate inequality method of `lhs`
            # (e.g., `lhs.__lt__`).  That method will try to reduce to
            # boolean or raise an exception.  It may keep calling
            # superclasses until it reaches `Expr` (e.g., `Expr.__lt__`).
            # In some cases, `Expr` will just invoke us again (if neither it
            # nor a subclass was able to reduce to boolean or raise an
            # exception).  In that case, it must call us with
            # `evaluate=False` to prevent infinite recursion.
            r = cls._eval_relation(lhs, rhs)
            if r is not None:
                return r
            # Note: not sure r could be None, perhaps we never take this
            # path?  In principle, could use this to shortcut out if a
            # class realizes the inequality cannot be evaluated further.

        # make a "non-evaluated" Expr for the inequality
        return Relational.__new__(cls, lhs, rhs, **options)


class _Greater(_Inequality):
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


class _Less(_Inequality):
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

    >>> e = x < y < z
    Traceback (most recent call last):
    ...
    TypeError: symbolic boolean expression has no truth value.

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
       http://docs.python.org/2/reference/expressions.html#notin).  This is done
       in an efficient way, so that each object being compared is only
       evaluated once and the comparison can short-circuit.  For example, ``1
       > 2 > 3`` is evaluated by Python as ``(1 > 2) and (2 > 3)``.  The
       ``and`` operator coerces each side into a bool, returning the object
       itself when it short-circuits.  The bool of the --Than operators
       will raise TypeError on purpose, because SymPy cannot determine the
       mathematical ordering of symbolic expressions.  Thus, if we were to
       compute ``x > y > z``, with ``x``, ``y``, and ``z`` being Symbols,
       Python converts the statement (roughly) into these steps:

        (1) x > y > z
        (2) (x > y) and (y > z)
        (3) (GreaterThanObject) and (y > z)
        (4) (GreaterThanObject.__nonzero__()) and (y > z)
        (5) TypeError

       Because of the "and" added at step 2, the statement gets turned into a
       weak ternary statement, and the first object's __nonzero__ method will
       raise TypeError.  Thus, creating a chained inequality is not possible.

           In Python, there is no way to override the ``and`` operator, or to
           control how it short circuits, so it is impossible to make something
           like ``x > y > z`` work.  There was a PEP to change this,
           :pep:`335`, but it was officially closed in March, 2012.

    .. [2] For more information, see these two bug reports:

       "Separate boolean and symbolic relationals"
       `Issue 4986 <https://github.com/sympy/sympy/issues/4986>`_

       "It right 0 < x < 1 ?"
       `Issue 6059 <https://github.com/sympy/sympy/issues/6059>`_

    """
    __slots__ = ()

    rel_op = '>='

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        # We don't use the op symbol here: workaround issue #7951
        return _sympify(lhs.__ge__(rhs))

Ge = GreaterThan


class LessThan(_Less):
    __doc__ = GreaterThan.__doc__
    __slots__ = ()

    rel_op = '<='

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        # We don't use the op symbol here: workaround issue #7951
        return _sympify(lhs.__le__(rhs))

Le = LessThan


class StrictGreaterThan(_Greater):
    __doc__ = GreaterThan.__doc__
    __slots__ = ()

    rel_op = '>'

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        # We don't use the op symbol here: workaround issue #7951
        return _sympify(lhs.__gt__(rhs))

Gt = StrictGreaterThan


class StrictLessThan(_Less):
    __doc__ = GreaterThan.__doc__
    __slots__ = ()

    rel_op = '<'

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        # We don't use the op symbol here: workaround issue #7951
        return _sympify(lhs.__lt__(rhs))

Lt = StrictLessThan


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
