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
