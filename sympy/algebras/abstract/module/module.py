from ..structure import AlgebraicStructure
from ..group import AbelianGroup
from ..ring import Ring

__all__ = [
    "Module",
]

class Module(AlgebraicStructure):
    """
    A base class for algebraic module.

    Explanation
    ===========

    Module is composite algebraic structure consists of a ring $R$ which is called scalar,
    an abelian group $M$, and scalar multiplication operation between $R$ and $M$.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of two structures
        The first one is scalar ring, and the
        second one is abelian group.

    operators : tuple of two Maps
        The first one is addition operator, and the
        second is multiplication operator.

    """
    def __new__(cls, name, sets, operators, **kwargs):
        if not len(sets) == 2:
            raise TypeError("%s must consist of two underlying structures." % cls)
        if not (len(operators) == 1 and all(o.arity == 2 for o in operators)):
            raise TypeError("%s must consist of one binary scalar multiplication operator." % cls)

        scalar_ring, abelian_group = sets

        if not isinstance(scalar_ring, Ring):
            raise TypeError("%s is not ring." % scalar_ring)
        if not isinstance(abelian_group, AbelianGroup):
            raise TypeError("%s is not abelian group." % abelian_group)

        return super().__new__(cls, name, sets, operators)

    @property
    def ring(self):
        return self.args[1].args[0]

    @property
    def group(self):
        return self.args[1].args[1]

    @property
    def scalar_multiplication(self):
        return self.args[2].args[0]
