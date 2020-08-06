from .ring import CommutativeRing
from ..group import AbelianGroup

__all__ = [
    "Field",
]

class Field(CommutativeRing):
    """
    A base class for algebraic field.

    Explanation
    ===========

    Field is a commutative ring whose domain and multiplication operator
    form Abelian group. Generalized division can be defined in ring.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets

    operators : tuple of two Maps
        The first one is addition operator, and the
        second is multiplication operator.

    Examples
    ========

    >>> from sympy import Field, Set, AdditionOperator, MultiplicationOperator

    Build a purely abstract field with no number involved.

    >>> A = Set('A')
    >>> a, b, e1, e2 = [A.element(n) for n in ('a', 'b', 'e1', 'e2')]

    >>> add = AdditionOperator(A**2, A, e1)
    >>> mul = MultiplicationOperator(A**2, A, e2)

    >>> F = Field('F', (A,), (add, mul))

    >>> F.div(a, b)
    a*(b**(-1))
    >>> F.div(a, a, evaluate=True)
    e2

    """
    def __new__(cls, name, sets, operators, **kwargs):
        obj = super().__new__(cls, name, sets, operators)

        add, mul = operators
        obj._add_group = AbelianGroup(name, sets, (add,))
        obj._mul_monoid = obj._mul_group = AbelianGroup(name, sets, (mul,))
        return obj

    @property
    def mul_group(self):
        return obj._mul_group

    @property
    def division_operator(self):
        return self.mul_op.right_division_operator()
    div_op = division_operator

    def divide(self, a, b, evaluate=False):
        return self.div_op(a, b, add_op=self.add_op, evaluate=evaluate)
    div = divide
