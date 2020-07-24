from sympy.assumptions import ask, Q
from .magma import Magma

__all__ = [
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
]

class LeftQuasigroup(Magma):
    r"""
    Left quasigroup is magma that left division operator can be derived from its operator.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, LeftQuasigroup

    >>> S = Set('S')
    >>> a, b, c = S.element('a'), S.element('b'), S.element('c')

    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = True
    >>> op = Op()

    >>> Q = LeftQuasigroup('Q', (S,), (op,))

    Left division exists.

    >>> op_ld = op.left_division()
    >>> op_ld(a, b)
    a \ b

    Divisions can be cancelled.

    >>> op(a, op_ld(a, b), evaluate=True)
    b
    >>> op_ld(a, op(a, b), evaluate=True)
    b

    """

    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.left_divisible(op)):
            raise TypeError("Left division of %s does not exist." % op)

        return super().__new__(cls, name, sets, operators)

class RightQuasigroup(Magma):
    r"""
    Right quasigroup is magma that right division operator can be derived from its operator.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, RightQuasigroup

    >>> S = Set('S')
    >>> a, b, c = S.element('a'), S.element('b'), S.element('c')

    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_right_divisible = True
    >>> op = Op()

    >>> Q = RightQuasigroup('Q', (S,), (op,))

    Left division exists.

    >>> op_rd = op.right_division()
    >>> op_rd(a, b)
    a / b

    Divisions can be cancelled.

    >>> op(op_rd(a, b), b, evaluate=True)
    a
    >>> op_rd(op(a, b), b, evaluate=True)
    a

    """

    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.right_divisible(op)):
            raise TypeError("Left division of %s does not exist." % op)

        return super().__new__(cls, name, sets, operators)

class Quasigroup(LeftQuasigroup, RightQuasigroup):
    r"""
    Quasigroup is both left quasigroup and right quasigroup.

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

    Left division and right division exist.

    >>> op_ld, op_rd = op.left_division(), op.right_division()
    >>> op_ld(a, b)
    a \ b
    >>> op_rd(a, b)
    a / b

    Divisions can be cancelled.

    >>> op(a, op_ld(a, b), evaluate=True)
    b
    >>> op_ld(a, op(a, b), evaluate=True)
    b
    >>> op(op_rd(a, b), b, evaluate=True)
    a
    >>> op_rd(op(a, b), b, evaluate=True)
    a

    """
