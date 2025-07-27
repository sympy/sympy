"""Implementation of :class:`CharacteristicZero` class. """


from sympy.polys.domains.domain import Domain
from sympy.utilities import public
from typing import Literal

@public
class CharacteristicZero(Domain):
    """Domain that has infinite number of elements. """

    has_CharacteristicZero = True

    def characteristic(self) -> Literal[0]:
        """Return the characteristic of this domain. """
        return 0
