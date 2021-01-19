from sympy.assumptions import Predicate
from sympy.multipledispatch import Dispatcher


class PrimePredicate(Predicate):
    """
    Prime number predicate.

    Explanation
    ===========

    ``ask(Q.prime(x))`` is true iff ``x`` is a natural number greater
    than 1 that has no positive divisors other than ``1`` and the
    number itself.

    Examples
    ========

    >>> from sympy import Q, ask
    >>> ask(Q.prime(0))
    False
    >>> ask(Q.prime(1))
    False
    >>> ask(Q.prime(2))
    True
    >>> ask(Q.prime(20))
    False
    >>> ask(Q.prime(-3))
    False

    """
    name = 'prime'
    handler = Dispatcher(
        "PrimeHandler",
        doc=("Handler for key 'prime'. Test that an expression represents a prime"
        " number. When the expression is an exact number, the result (when True)"
        " is subject to the limitations of isprime() which is used to return the "
        "result.")
    )
