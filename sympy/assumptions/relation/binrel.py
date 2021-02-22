"""
General binary relations.
"""

from sympy import S
from sympy.assumptions import AppliedPredicate, ask, Predicate
from sympy.core.compatibility import ordered
from sympy.core.kind import BooleanKind
from sympy.core.parameters import global_parameters
from sympy.logic.boolalg import BooleanAtom

__all__ = ["BinaryRelation", "AppliedBinaryRelation"]


class BinaryRelation(Predicate):
    """
    Base class for all binary relational predicates.

    Explanation
    ===========

    Binary relation takes two arguments and returns ``AppliedBinaryRelation``
    instance. To evaluate it to boolean value, use :obj:`~.ask()` or
    :obj:`~.refine()` function.

    You can add support for new types by registering the handler to dispatcher.
    See :obj:`~.Predicate()` for more information about predicate dispatching.

    Examples
    ========

    Applying and evaluating to boolean value:

    >>> from sympy import Q, ask, sin, cos
    >>> from sympy.abc import x
    >>> Q.eq(sin(x)**2+cos(x)**2, 1)
    sin(x)**2 + cos(x)**2 = 1
    >>> ask(_)
    True

    You can define a new binary relation by subclassing and dispatching.
    Here, we define a relation $R$ such that $x R y$ returns true if
    $x = y + 1$.

    >>> from sympy import ask, Number
    >>> from sympy.assumptions import BinaryRelation
    >>> class MyRel(BinaryRelation):
    ...     name = "R"
    ...     is_reflexive = False
    >>> R = MyRel()
    >>> @R.register(Number, Number)
    ... def _(n1, n2, assumptions):
    ...     return ask(Q.zero(n1 - n2 - 1), assumptions)
    >>> R(2, 1)
    2 R 1

    Now, we can use ``ask()`` to evaluate it to boolean value.

    >>> ask(R(2, 1))
    True
    >>> ask(R(1, 2))
    False

    ``R`` returns ``False`` with minimum cost if two arguments have same
    structure because R is antireflexive relation [1] by
    ``is_reflexive = False``.

    >>> ask(R(x, x))
    False

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Reflexive_relation
    """

    is_reflexive = None
    is_symmetric = None

    is_Relational = True    # compatibility for old Relational

    def __call__(self, *args):
        if not len(args) == 2:
            raise ValueError("Binary relation takes two arguments, but got %s." % len(args))
        return AppliedBinaryRelation(self, *args)

    @property
    def reversed(self):
        if self.is_symmetric:
            return self
        return None

    @property
    def negated(self):
        return None

    def _compare_reflexive(self, lhs, rhs):
        # quick exit for structurally same arguments
        # do not check != here because it cannot catch the
        # equivalent arguements with different structures.

        # reflexivity does not hold to NaN
        if lhs is S.NaN or rhs is S.NaN:
            return None

        reflexive = self.is_reflexive
        if reflexive is None:
            pass
        elif reflexive and (lhs == rhs):
            return True
        elif not reflexive and (lhs == rhs):
            return False
        return None

    def eval(self, args, assumptions=True):
        # quick exit for structurally same arguments
        ret = self._compare_reflexive(*args)
        if ret is not None:
            return ret

        # don't perform simplify here. (done by AppliedBinaryRelation._eval_ask)
        # evaluate by multipledispatch
        lhs, rhs = args
        ret = self.handler(lhs, rhs, assumptions=assumptions)
        if ret is not None:
            return ret

        # check reversed order if the relation is reflexive
        if self.is_reflexive:
            types = (type(lhs), type(rhs))
            if self.handler.dispatch(*types) is not self.handler.dispatch(*reversed(types)):
                ret = self.handler(rhs, lhs, assumptions=assumptions)

        return ret

    def _simplify_applied(self, lhs, rhs, **kwargs):
        lhs, rhs = lhs.simplify(**kwargs), rhs.simplify(**kwargs)
        return self(lhs, rhs)

    def _eval_binary_symbols(self, lhs, rhs):
        # override where necessary
        return set()


class AppliedBinaryRelation(AppliedPredicate):
    """
    The class of expressions resulting from applying ``BinaryRelation``
    to the arguments.

    """

    @property
    def lhs(self):
        """The left-hand side of the relation."""
        return self.arguments[0]

    @property
    def rhs(self):
        """The right-hand side of the relation."""
        return self.arguments[1]

    @property
    def reversed(self):
        """
        Try to return the relationship with sides reversed.
        """
        revfunc = self.function.reversed
        if revfunc is None:
            return self
        return revfunc(self.rhs, self.lhs)

    @property
    def reversedsign(self):
        """
        Try to return the relationship with signs reversed.
        """
        revfunc = self.function.reversed
        if revfunc is None:
            return self
        if not any(side.kind is BooleanKind for side in self.arguments):
            return revfunc(-self.lhs, -self.rhs)
        return self

    @property
    def negated(self):
        neg_rel = self.function.negated
        if neg_rel is None:
            return ~self
        return neg_rel(*self.arguments)

    @property
    def canonical(self):
        """
        Return a canonical form of the relational by putting a
        number on the rhs, canonically removing a sign or else
        ordering the args canonically. No other simplification is
        attempted.
        """
        args = self.arguments
        r = self
        if r.rhs.is_number:
            if any(side is S.NaN for side in (r.lhs, r.rhs)):
                pass
            elif r.rhs.is_Number and r.lhs.is_Number and r.lhs > r.rhs:
                r = r.reversed
        elif r.lhs.is_number:
            r = r.reversed
        elif tuple(ordered(args)) != args:
            r = r.reversed

        LHS_CEMS = getattr(r.lhs, 'could_extract_minus_sign', None)
        RHS_CEMS = getattr(r.rhs, 'could_extract_minus_sign', None)

        if any(side.kind is BooleanKind for side in r.arguments):
            return r

        # Check if first value has negative sign
        if LHS_CEMS and LHS_CEMS():
            return r.reversedsign
        elif not r.rhs.is_number and RHS_CEMS and RHS_CEMS():
            # Right hand side has a minus, but not lhs.
            # How does the expression with reversed signs behave?
            # This is so that expressions of the type
            # Eq(x, -y) and Eq(-x, y)
            # have the same canonical representation
            expr1, _ = ordered([r.lhs, -r.rhs])
            if expr1 != r.lhs:
                return r.reversed.reversedsign

        return r

    def _eval_ask(self, assumptions):
        ret = self.function.eval(self.arguments, assumptions)
        if ret is not None:
            return ret

        # simplify and try again
        rel = self.simplify()
        return rel.function.eval(rel.arguments, assumptions)

    def _eval_simplify(self, **kwargs):
        return self.function._simplify_applied(self.lhs, self.rhs, **kwargs)

    def __bool__(self):
        ret = ask(self)
        if ret is None:
            raise TypeError("Cannot determine truth value of %s" % self)
        return ret

    @property
    def binary_symbols(self):
        return self.function._eval_binary_symbols(*self.arguments)


class _DeprecatedRelational(AppliedBinaryRelation):
    """
    Class to make migration from ``core/relational`` to this module.
    Old signatures such as ``Eq(x, y)`` returns this class. Evaluates by
    default and by simplification.

    """
    def __new__(cls, predicate, *args, **options):
        evaluate = options.pop('evaluate', global_parameters.evaluate)
        if evaluate:
            val = predicate.eval(args)
            if val is not None:
                return val
        return super().__new__(cls, predicate, *args)

    def _eval_simplify(self, **kwargs):
        rel = super()._eval_simplify(**kwargs)
        ret = rel.refine()
        if isinstance(ret, BooleanAtom):
            return ret
        return rel
