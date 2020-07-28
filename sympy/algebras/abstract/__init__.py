__all__ = [
    "AlgebraicStructure",
    "Magma",
    "Semigroup",
    "LeftQuasigroup", "RightQuasigroup", "Quasigroup",
    "Monoid",
    "Loop",
    "Group", 'AbelianGroup',
    "Ring", "CommutativeRing",
    "Field",
    "Module",
    "VectorSpace",
]

from .structure import AlgebraicStructure

from .group import (
    Magma, Semigroup,
    LeftQuasigroup, RightQuasigroup, Quasigroup, Monoid, Loop,
    Group, AbelianGroup,
)

from .ring import (
    Ring, CommutativeRing, Field,
)

from .module import (
    Module, VectorSpace,
)

### Defined structures

from sympy.core.singleton import S
from sympy.map.add import AdditionOperator
from sympy.map.mul import MultiplicationOperator

int_add = AdditionOperator(S.Integers**2, S.Integers, S.Zero)
int_mul = MultiplicationOperator(S.Integers**2, S.Integers, S.One)
S.Integers_ring = CommutativeRing('Z', (S.Integers,), (int_add, int_mul))
S.Integers_ring._latex = lambda self, *args, **kwargs: r'\mathbb{Z}'

real_add = AdditionOperator(S.Reals**2, S.Reals, S.Zero)
real_mul = MultiplicationOperator(S.Reals**2, S.Reals, S.One)
S.Reals_field = Field('R', (S.Reals,), (real_add, real_mul))
S.Reals_field._latex = lambda self, *args, **kwargs: r'\mathbb{R}'

complexes_add = AdditionOperator(S.Complexes**2, S.Complexes, S.Zero)
complexes_mul = MultiplicationOperator(S.Complexes**2, S.Complexes, S.One)
S.Complexes_field = Field('C', (S.Complexes,), (complexes_add, complexes_mul))
S.Complexes_field._latex = lambda self, *args, **kwargs: r'\mathbb{C}'
