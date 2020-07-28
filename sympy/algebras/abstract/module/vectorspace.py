from ..ring import Field
from .module import Module

__all__ = [
    "VectorSpace",
]

class VectorSpace(Module):
    """
    A base class for vector space.

    Explanation
    ===========

    Vector space is a module whose ring is field.

    """
    def __new__(cls, name, sets, operators, **kwargs):
        scalar_field, abelian_group = sets

        if not isinstance(scalar_field, Field):
            raise TypeError("%s is not field." % scalar_field)

        return super().__new__(cls, name, sets, operators)

    @property
    def field(self):
        return self.args[1].args[0]
    scalars = field

    @property
    def vectors(self):
        return self.args[1].args[1]
