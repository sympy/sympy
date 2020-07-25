from sympy.assumptions import ask, Q
from .semigroup import Semigroup

__all__ = [
    "Monoid",
]

class Monoid(Semigroup):
    """
    Monoid is semigroup whose operation has two-sided identity.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Monoid

    >>> S = Set('S')
    >>> a, b, c, e = [S.element(n) for n in 'abce']

    >>> class MonoidOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_associative = True
    ...     identity = e
    >>> op = MonoidOp()

    >>> M = Monoid('M', (S,), (op,))

    Operation of monoid is associative.

    >>> op(a, op(b, c), evaluate=True)
    a * b * c

    Operation of monoid has identity.

    >>> op(a, op(e, b), evaluate=True)
    a * b

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if getattr(op, 'identity', None) is None:
            raise TypeError("%s does not have identity." % op)

        return super().__new__(cls, name, sets, operators)
