"""Abstract :class:`RingExtension` domain type. """

from __future__ import annotations

from abc import abstractmethod
from typing import Generic
from sympy.polys.domains.domain import Domain, Er, Eg

from sympy.utilities import public


@public
class RingExtension(Domain[Er], Generic[Er, Eg]):
    """Represents a ring extension.

    Examples are polynomial rings :ref:`K[x]`, algebraic number fields
    :ref:`QQ(a)`, the Guassian domains :ref:`ZZ_I` and :ref:`QQ_I`.

    A ``RingExtension`` domain has a ground domain and a set of generators
    and provides a ``.to_dict()`` method to convert elements to a dict
    representation over the ground domain.

    Examples
    ========

    >>> from sympy import QQ, Symbol
    >>> x = Symbol('x')
    >>> R = QQ[x]
    >>> R.to_dict(R(2))
    {(0,): 2}

    See Also
    ========

    AlgebraicField
    PolynomialRing
    GaussianIntegerRing
    GaussianRationalField
    """

    is_RingExtension = True

    dom: Domain[Eg]
    """The ground domain of the ring extension."""

    ngens: int
    """The number of generators of the ring extension."""

    # XXX: For polynomial rings gens is a tuple of domain elements whereas for
    # algebraic number fields gens is a tuple of Expr.
    # gens: tuple[ER, ...]

    @abstractmethod
    def to_dict(self, element: Er) -> dict[tuple[int, ...], Eg]:
        """Convert ``element`` to a dict over the ground domain."""
        raise NotImplementedError
