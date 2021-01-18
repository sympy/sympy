from sympy.assumptions import Predicate
from sympy.multipledispatch import Dispatcher

class NegativePredicate(Predicate):
    r"""
    Negative number predicate.

    Explanation
    ===========

    ``Q.negative(x)`` is true iff ``x`` is a real number and :math:`x < 0`, that is,
    it is in the interval :math:`(-\infty, 0)`.  Note in particular that negative
    infinity is not negative.

    A few important facts about negative numbers:

    - Note that ``Q.nonnegative`` and ``~Q.negative`` are *not* the same
        thing. ``~Q.negative(x)`` simply means that ``x`` is not negative,
        whereas ``Q.nonnegative(x)`` means that ``x`` is real and not
        negative, i.e., ``Q.nonnegative(x)`` is logically equivalent to
        ``Q.zero(x) | Q.positive(x)``.  So for example, ``~Q.negative(I)`` is
        true, whereas ``Q.nonnegative(I)`` is false.

    - See the documentation of ``Q.real`` for more information about
        related facts.

    Examples
    ========

    >>> from sympy import Q, ask, symbols, I
    >>> x = symbols('x')
    >>> ask(Q.negative(x), Q.real(x) & ~Q.positive(x) & ~Q.zero(x))
    True
    >>> ask(Q.negative(-1))
    True
    >>> ask(Q.nonnegative(I))
    False
    >>> ask(~Q.negative(I))
    True

    """
    name = 'negative'
    handler = Dispatcher(
        "NegativeHandler",
        doc=("Handler for Q.negative. Test that an expression is strictly less"
        " than zero.")
    )
