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
    form Abelian group.

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
        return self.div_op(a, b, evaluate=evaluate)
    div = divide

    @property
    def power_operator(self):
        return self.mul_op.exponent_operator()
    pow_op = power_operator

    def power(self, x, n, evaluate=False):
        return self.pow_op(x, n, evaluate=evaluate)
    pow = power
