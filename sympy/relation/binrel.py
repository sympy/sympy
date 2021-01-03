"""
Module to implement general binary relations.
"""

from sympy.assumptions import AppliedPredicate, Predicate
from sympy.core.compatibility import ordered
from sympy.logic.boolalg import BooleanAtom
from sympy.simplify import simplify


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
    >>> ask(Q.eq(sin(x)**2+cos(x)**2, 1))
    True

    You can define a new binary relation by subclassing and dispatching.
    Here, we define a relation $R$ such that $x R y$ returns true if
    $x = y + 1$.

    >>> from sympy import ask, Number
    >>> from sympy.relation import BinaryRelation
    >>> class MyRel(BinaryRelation):
    ...     name = "R"
    ...     is_reflexive = False
    >>> R = MyRel()
    >>> @R.register(Number, Number)
    ... def _(n1, n2, assumptions):
    ...     return ask(Q.zero(n1 - n2 - 1), assumptions)
    >>> R(2, 1)
    2 R 1

    Now, we can use ``ask()`` to evaluate the applied to boolean value.

    >>> ask(R(2, 1))
    True
    >>> ask(R(1, 2))
    False

    Although we didn't dispatch the types other than ``Number``, ``R``
    returns ``False`` with minimum cost if two arguments have same tree
    structure because R is antireflexive relation [1].

    >>> ask(R(x, x))
    False

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Reflexive_relation
    """

    is_reflexive = None
    is_symmetric = None

    def __call__(self, *args):
        if not len(args) == 2:
            raise ValueError("Binary relation takes two arguments, but got %s." % len(args))
        return AppliedBinaryRelation(self, *args)

    @property
    def reversed(self):
        return None

    def _compare_reflexive(self, lhs, rhs):
        # quick exit for structurally same arguments
        # do not check != here because it cannot catch the
        # equivalent arguements with different structures.
        reflexive = self.is_reflexive
        if reflexive is None:
            pass
        elif reflexive and (lhs == rhs):
            return True
        elif not reflexive and (lhs == rhs):
            return False
        return None


class AppliedBinaryRelation(AppliedPredicate):
    """
    The class of expressions resulting from applying ``BinaryRelation``
    to the arguments.

    This class wraps its argument and remain unevaluated. To evaluate it
    to boolean value, use :obj:`~.ask()` function.

    Examples
    ========

    >>> from sympy import cos, sin, Q
    >>> from sympy.abc import x, y, z
    >>> eqn1 = Q.eq(sin(x)**2 + cos(x)**2, 1)
    >>> eqn1
    sin(x)**2 + cos(x)**2 = 1

    ``.simplify()`` simplifies the relation, but does not evaluate it
    even if the relation is identical. Also, it does not take assumption
    into account.

    >>> eqn1.simplify()
    0 = 0
    >>> eqn2 = Q.eq(x*y, x*z)
    >>> eqn2.simplify()
    x*y = x*z

    ``ask()`` evaluates the relation to boolean value. If the truth
    value cannot be determined, it returns ``None``.

    >>> from sympy import ask, Abs
    >>> eqn3 = Q.eq(Abs(x), x)
    >>> print(ask(eqn2))
    None
    >>> ask(eqn3, Q.positive(x))
    True
    >>> ask(eqn3, Q.negative(x))
    False

    """

    # will be deleted after _op_priority is removed from SymPy
    _op_priority = 1000

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
        """Return the relationship with sides reversed.

        """
        revfunc = self.function.reversed
        if revfunc is None:
            if self.function.is_symmetric:
                return self.function(self.lhs, self.rhs)
            return self
        return self.function.reversed(self.rhs, self.lhs)

    @property
    def reversedsign(self):
        """Return the relationship with signs reversed.

        """
        a, b = self.arguments
        if not (isinstance(a, BooleanAtom) or isinstance(b, BooleanAtom)):
            return self.function.reversed(-self.lhs, -self.rhs)
        else:
            return self

    @property
    def canonical(self):
        """Return a canonical form of the relational by putting a
        number on the rhs, canonically removing a sign or else
        ordering the args canonically. No other simplification is
        attempted.

        """
        args = self.arguments
        r = self
        if r.rhs.is_number:
            if r.rhs.is_Number and r.lhs.is_Number and r.lhs > r.rhs:
                r = r.reversed
        elif r.lhs.is_number:
            r = r.reversed
        elif tuple(ordered(args)) != args:
            r = r.reversed

        LHS_CEMS = getattr(r.lhs, 'could_extract_minus_sign', None)
        RHS_CEMS = getattr(r.rhs, 'could_extract_minus_sign', None)

        if isinstance(r.lhs, BooleanAtom) or isinstance(r.rhs, BooleanAtom):
            return r

        # Check if first value has negative sign
        if LHS_CEMS and LHS_CEMS():
            return r.reversedsign
        elif not r.rhs.is_number and RHS_CEMS and RHS_CEMS():
            # Right hand side has a minus, but not lhs.
            # How does the expression with reversed signs behave?
            # This is so that expressions of the type
            # Q.eq(x, -y) and Q.eq(-x, y)
            # have the same canonical representation
            expr1, _ = ordered([r.lhs, -r.rhs])
            if expr1 != r.lhs:
                return r.reversed.reversedsign

        return r

    @property
    def binary_symbols(self):
        return set()

    def as_Relational(self):
        return self.function.as_Relational(*self.arguments)

    def _eval_ask(self, assumptions):
        # quick exit for structurally same arguments
        ret = self.function._compare_reflexive(self.lhs, self.rhs)
        if ret is not None:
            return ret

        lhs, rhs = self.lhs.refine(assumptions), self.rhs.refine(assumptions)
        r = self.function(lhs, rhs).simplify()

        # attempt quick exit again with simplified arguments
        ret = r.function._compare_reflexive(r.lhs, r.rhs)
        if ret is not None:
            return ret

        # use multipledispatch handler
        return r.function.eval(r.arguments, assumptions)

    def simplify(self, equation=True, side="all", **kwargs):
        """
        Simplify *self* without evaluating to boolean value. Assumption
        is not taken into consideration.

        Parameters
        ==========

        equation : bool, optional
            If ``True``, simplify both sides and canonical result.
            If ``False``, each side is simplified but not canonicalized.
            You can decide which side to be simplified by *side* argument.

        side : "all", "lhs" or "rhs", optional
            Specify which side the rewriting will be done. Only valid
            when *equation* is ``True``. Default is "all".

        Examples
        ========

        >>> from sympy import cos, gamma, sin, Q
        >>> from sympy.abc import x
        >>> eqn = Q.eq(sin(x)**2 + cos(x)**2, gamma(x)/gamma(x-2))
        >>> eqn.simplify()
        x**2 - 3*x = -1
        >>> eqn.simplify(equation=False)
        1 = (x - 2)*(x - 1)
        >>> eqn.simplify(equation=False, side='lhs')
        1 = gamma(x)/gamma(x - 2)
        """
        kwargs.update(equation=equation,
                      side=side)
        return simplify(self, **kwargs)

    def _eval_simplify(self, **kwargs):
        from .eqntools import eqnsimp

        equation = kwargs.get('equation', True)
        side = kwargs.get('side', 'all')

        lhs, rhs = self.arguments

        if equation:
            return eqnsimp(self.function, lhs, rhs, **kwargs)

        if side in ("all", "lhs"):
            lhs = self.lhs.simplify(**kwargs)
        if side in ("all", "rhs"):
            rhs = self.rhs.simplify(**kwargs)
        return self.function(lhs, rhs)
