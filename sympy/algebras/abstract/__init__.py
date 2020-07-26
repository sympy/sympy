__all__ = [
    "AlgebraicStructure",
    "Magma",
    "Semigroup",
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
    "Monoid",
    "Loop",
    "Group", 'AbelianGroup',
]

from .structure import AlgebraicStructure

from .group import (
    Magma, Semigroup,
    LeftQuasigroup, RightQuasigroup, Quasigroup, Monoid, Loop,
    Group, AbelianGroup,
)
