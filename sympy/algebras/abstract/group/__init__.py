"""Group-like structure"""

__all__ = [
    "Magma",
    "Semigroup",
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
    "Monoid",
    "Loop",
    "Group", 'AbelianGroup',
]

from .magma import Magma

from .semigroup import Semigroup

from .quasigroup import LeftQuasigroup, RightQuasigroup, Quasigroup

from .monoid import Monoid

from .loop import Loop

from .group import Group, AbelianGroup
