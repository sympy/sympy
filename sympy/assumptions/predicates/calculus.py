from sympy.assumptions import Predicate
from sympy.multipledispatch import Dispatcher

class FinitePredicate(Predicate):
    """
    Finite predicate.

    Explanation
    ===========

    ``Q.finite(x)`` is true if ``x`` is neither an infinity
    nor a ``NaN``. In other words, ``ask(Q.finite(x))`` is true for all ``x``
    having a bounded absolute value.

    Examples
    ========

    >>> from sympy import Q, ask, Symbol, S, oo, I
    >>> x = Symbol('x')
    >>> ask(Q.finite(S.NaN))
    False
    >>> ask(Q.finite(oo))
    False
    >>> ask(Q.finite(1))
    True
    >>> ask(Q.finite(2 + 3*I))
    True
    >>> print(ask(Q.finite(x), Q.positive(x)))
    None

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Finite

    """
    name = 'finite'
    handler = Dispatcher(
        "FiniteHandler",
        doc=("Handler for Q.finite. Test that an expression is bounded respect"
        " to all its variables.")
    )


class InfinitePredicate(Predicate):
    """
    Infinite number predicate.

    ``Q.infinite(x)`` is true iff the absolute value of ``x`` is
    infinity.

    """
    # TODO: Add examples
    name = 'infinite'
    handler = Dispatcher(
        "InfiniteHandler",
        doc="""Handler for Q.infinite key."""
    )
