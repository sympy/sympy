"""Implementaton of :class:`CharacteristicZero` class. """

__all__ = ["CharacteristicZero"]

from sympy.polys.domains.domain import Domain


class CharacteristicZero(Domain):
    """Domain that has infinite number of elements. """

    has_CharacteristicZero = True

    def characteristic(self):
        """Return the characteristic of this domain. """
        return 0
