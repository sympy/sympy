from sympy.assumptions import ask, Q
from .magma import Magma

__all__ = [
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
]

class LeftQuasigroup(Magma):
    r"""
    A base class for left quasigroup.

    Explanation
    ===========

    Left quasigroup is magma that left division operator can be derived from its operator.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets
        See sympy.sets module.

    operators : tuple of one Map
        See sympy.map module.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, LeftQuasigroup

    >>> S = Set('S')
    >>> a, b, c = [S.element(n) for n in 'abc']

    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = True
    >>> op = Op()

    >>> Q = LeftQuasigroup('Q', (S,), (op,))
    >>> Q_op = Q.operator

    Left division exists.

    >>> Q_ld = Q.left_division
    >>> Q_ld(a, b)
    a \ b

    Divisions can be cancelled.

    >>> Q_op(a, Q_ld(a, b), evaluate=True)
    b
    >>> Q_ld(a, Q_op(a, b), evaluate=True)
    b

    """

    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.left_divisible(op)):
            raise TypeError("Left division of %s does not exist." % op)

        return super().__new__(cls, name, sets, operators)

    @property
    def left_division(self):
        return self.operator.left_division_operator()

class RightQuasigroup(Magma):
    r"""
    A base class for righ quasigroup.

    Explanation
    ===========

    Right quasigroup is magma that right division operator can be derived from its operator.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets
        See sympy.sets module.

    operators : tuple of one Map
        See sympy.map module.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, RightQuasigroup

    >>> S = Set('S')
    >>> a, b, c = [S.element(n) for n in 'abc']

    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_right_divisible = True
    >>> op = Op()

    >>> Q = RightQuasigroup('Q', (S,), (op,))
    >>> Q_op = Q.operator

    Left division exists.

    >>> Q_rd = Q.right_division
    >>> Q_rd(a, b)
    a / b

    Divisions can be cancelled.

    >>> Q_op(Q_rd(a, b), b, evaluate=True)
    a
    >>> Q_rd(Q_op(a, b), b, evaluate=True)
    a

    """

    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.right_divisible(op)):
            raise TypeError("Left division of %s does not exist." % op)

        return super().__new__(cls, name, sets, operators)

    @property
    def right_division(self):
        return self.operator.right_division_operator()

class Quasigroup(LeftQuasigroup, RightQuasigroup):
    r"""
    A base class for quasigroup.

    Explanation
    ===========

    Quasigroup is both left quasigroup and right quasigroup.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets
        See sympy.sets module.

    operators : tuple of one Map
        See sympy.map module.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Quasigroup

    >>> S = Set('S')
    >>> a, b, c = [S.element(n) for n in 'abc']

    >>> class QuasigroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = is_right_divisible = True
    >>> op = QuasigroupOp()

    >>> Q = Quasigroup('Q', (S,), (op,))
    >>> Q_op = Q.operator

    Left division and right division exist.

    >>> Q_ld, Q_rd = Q.left_division, Q.right_division
    >>> Q_ld(a, b)
    a \ b
    >>> Q_rd(a, b)
    a / b

    Divisions can be cancelled.

    >>> Q_op(a, Q_ld(a, b), evaluate=True)
    b
    >>> Q_ld(a, Q_op(a, b), evaluate=True)
    b
    >>> Q_op(Q_rd(a, b), b, evaluate=True)
    a
    >>> Q_rd(Q_op(a, b), b, evaluate=True)
    a

    """
    # All attributes are inherited
    pass
