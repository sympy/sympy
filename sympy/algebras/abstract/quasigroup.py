from sympy.assumptions import ask, Q
from .magma import Magma

__all__ = [
    "Quasigroup",
]

class Quasigroup(Magma):
    r"""
    Quasigroup is magma whose operation can be defined with left and right division.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Quasigroup

    >>> S = Set('S')
    >>> a, b, c = S.element('a'), S.element('b'), S.element('c')

    >>> class QuasigroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = is_right_divisible = True
    >>> op = QuasigroupOp()

    >>> Q = Quasigroup('Q', (S,), (op,))

    >>> op(a, b)
    a * b
    >>> op(a, b) in Q
    True

    Left division and right division exists

    >>> op_ld, op_rd = op.left_division(), op.right_division()
    >>> op_ld(a, b)
    a \ b
    >>> op_rd(a, b)
    a / b

    Divisions can be cancelled

    >>> op(a, op_ld(a, b), evaluate=True)
    b
    >>> op_ld(a, op(a, b), evaluate=True)
    b
    >>> op(op_rd(a, b), b, evaluate=True)
    a
    >>> op_rd(op(a, b), b, evaluate=True)
    a

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.left_divisible(op)):
            raise TypeError("Left division of %s does not exist." % op)
        if not ask(Q.right_divisible(op)):
            raise TypeError("Left division of %s does not exist." % op)

        return super().__new__(cls, name, sets, operators)
