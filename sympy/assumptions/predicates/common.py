from sympy.assumptions import Predicate
from sympy.multipledispatch import Dispatcher


class CommutativePredicate(Predicate):
    """
    Commutative predicate.

    Explanation
    ===========

    ``ask(Q.commutative(x))`` is true iff ``x`` commutes with any other
    object with respect to multiplication operation.

    """
    # TODO: Add examples
    name = 'commutative'
    handler = Dispatcher("FiniteHandler", doc="Handler for key 'commutative'.")


class IsTruePredicate(Predicate):
    """
    Generic predicate.

    Explanation
    ===========

    ``ask(Q.is_true(x))`` is true iff ``x`` is true. This only makes
    sense if ``x`` is a predicate.

    Examples
    ========

    >>> from sympy import ask, Q, symbols
    >>> x = symbols('x')
    >>> ask(Q.is_true(True))
    True

    """
    name = 'is_true'
    handler = Dispatcher(
        "FiniteHandler",
        doc="Wrapper allowing to query the truth value of a boolean expression."
    )

    def __call__(self, *args):
        if len(args) == 1 and hasattr(args[0], "as_Predicate"):
            return args[0].as_Predicate()
        return super().__call__(*args)
