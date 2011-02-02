from sympy import Expr

from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.commutator import Commutator

def test_commutator_dagger():
    A = Expr('A', **{'commutative':False})
    B = Expr('B', **{'commutative':False})
    C = Expr('C', **{'commutative':False})
    comm = Commutator(A*B,C)
    assert Dagger(comm).expand(commutator=True) ==\
        - Commutator(Dagger(B),Dagger(C))*Dagger(A) -\
        Dagger(B)*Commutator(Dagger(A),Dagger(C))