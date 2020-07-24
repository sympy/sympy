__all__ = [
    "AlgebraicStructure",
    "Magma",
    "Semigroup",
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
    "Monoid",
]

from .structure import AlgebraicStructure

from .magma import Magma

from .semigroup import Semigroup

from .quasigroup import LeftQuasigroup, RightQuasigroup, Quasigroup

from .monoid import Monoid
