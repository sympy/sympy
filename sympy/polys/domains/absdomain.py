"""Abstract domain supporting absolute values. """

from __future__ import annotations

from abc import abstractmethod
from typing import Generic

from sympy.polys.domains.domain import Domain, Er, Eordered
from sympy.utilities import public


@public
class AbsDomain(Domain[Er], Generic[Er, Eordered]):
    """A domain whose absolute values belong to an ordered domain.

    The element type and the absolute value type can differ. For example,
    absolute values in :ref:`CC` belong to :ref:`RR` at the same precision.
    For :ref:`ZZ`, :ref:`QQ` and :ref:`RR`, absolute values belong to the
    original domain.

    Examples
    ========

    >>> from sympy import CC, ZZ
    >>> ZZ.abs_domain == ZZ
    True
    >>> CC.abs_domain
    RR
    >>> CC.abs(CC(3, 4)) == CC.abs_domain(5)
    True
    """

    is_AbsDomain = True

    @property
    @abstractmethod
    def abs_domain(self) -> Domain[Eordered]:
        """The ordered domain containing the absolute values."""
        raise NotImplementedError

    @abstractmethod
    def abs(self, a: Er) -> Eordered:
        """Return the absolute value of ``a`` in :attr:`abs_domain`."""
        raise NotImplementedError
