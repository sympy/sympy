from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.state import Bra, Ket

def test_innerproduct_dagger():
    k = Ket('k')
    b = Bra('b')
    ip = b*k
    assert Dagger(ip) == Dagger(k)*Dagger(b)