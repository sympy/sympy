from sympy.assumptions import ask, Q
from .magma import Magma

__all__ = [
    "Semigroup",
]

class Semigroup(Magma):
    """
    Semigroup is magma whose operation is associative.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Semigroup

    >>> S = Set('S')
    >>> a, b, c = S.element('a'), S.element('b'), S.element('c')

    >>> class SemigroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_associative = True
    >>> op = SemigroupOp()

    >>> G = Semigroup('G', (S,), (op,))

    Operation of semigroup is associative.

    >>> op(a, op(b, c), evaluate=True)
    a * b * c

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.associative(op)):
            raise TypeError("%s is not associative." % op)

        return super().__new__(cls, name, sets, operators)
