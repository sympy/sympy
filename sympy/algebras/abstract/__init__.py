__all__ = [
    "AlgebraicStructure",
    "Magma",
    "Semigroup",
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
    "Monoid",
    "Loop",
    "Group", 'AbelianGroup',
    "Ring", "CommutativeRing",
]

from .structure import AlgebraicStructure

from .group import (
    Magma, Semigroup,
    LeftQuasigroup, RightQuasigroup, Quasigroup, Monoid, Loop,
    Group, AbelianGroup,
)

from .ring import (
    Ring, CommutativeRing,
)
