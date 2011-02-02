from sympy import Expr

from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.anticommutator import AntiCommutator


def test_anticommutator_dagger():
    A = Expr('A', **{'commutative':False})
    B = Expr('C', **{'commutative':False})
    assert Dagger(AntiCommutator(A, B)) == AntiCommutator(Dagger(A),Dagger(B))