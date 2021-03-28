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
    handler = Dispatcher("CommutativeHandler", doc="Handler for key 'commutative'.")


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
        "IsTrueHandler",
        doc="Wrapper allowing to query the truth value of a boolean expression."
    )
