from __future__ import annotations

from .cartan_type import CartanType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.matrices import Matrix

def CartanMatrix(ct: str | list[str | int]) -> Matrix:
    """Access the Cartan matrix of a specific Lie algebra

    Examples
    ========

    >>> from sympy.liealgebras.cartan_matrix import CartanMatrix
    >>> CartanMatrix("A2")
    Matrix([
    [ 2, -1],
    [-1,  2]])

    >>> CartanMatrix(['C', 3])
    Matrix([
    [ 2, -1,  0],
    [-1,  2, -1],
    [ 0, -2,  2]])

    This method works by returning the Cartan matrix
    which corresponds to Cartan type t.
    """

    return CartanType(ct).cartan_matrix()
