from sympy.physics.quantum import (
    Operator,
    Dagger,
    State,
    Ket,
    Bra,
    TimeDepKet,
    TimeDepBra,
    BasisSet,
    InnerProduct,
    OuterProduct
)
from sympy.physics.hilbert import (
    HilbertSpace,
    l2, 
    L2, 
    FockSpace, 
    TensorProductHilbertSpace, 
    DirectSumHilbertSpace, 
    TensorPowerHilbertSpace
)
from sympy.physics.quantumbasic import QuantumError, QuantumBasic
from sympy import Symbol, expand
from sympy.utilities.pytest import XFAIL

def test_QMul_Ket():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    assert a*A
    assert a*A == A*a
    assert expand(a*(A+B)) == a*A+a*B

@XFAIL
def test_KetKet():
    A = Ket('a')
    B = Bra('b')
    assert A*B

def test_QMul_Bra():
    a = Symbol('a')
    b = Symbol('b')
    C = Bra('C')
    D = Bra('D')
    assert a*C
    assert a*C == C*a
    assert expand(a*(C+D)) == a*C+a*D

def test_QMul_Operator():
    a = Symbol('a')
    b = Symbol('b')
    E = Operator('E')
    F = Operator('F')
    G = Operator('G')
    assert a*E
    assert expand(a*(E+F)) == a*E+a*F
    assert E*F != F*E
    assert E*(F*G) == (E*F)*G
    assert Dagger(a*E) == Dagger(E)*Dagger(a)

def test_QMul_mixed():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    E = Operator('E')
    F = Operator('F')
    G = Operator('G')
    assert expand((E+F)*(A+B)) == (E*A + E*B + F*A + F*B)
    assert E*B
    assert Dagger(E*A) == Dagger(A)*Dagger(E) == Bra('A')*Dagger(C)

def test_QAdd():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    C = Bra('C')
    D = Bra('D')
    E = Operator('E')
    F = Operator('F')
    G = Operator('G')
    assert A+B
    assert A+B == B+A
    assert C+D
    assert C+D == D+C
    assert E+F == F+E
    assert (E+F)+G == E+(F+G)

def test_QPow():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    E = Operator('E')
    F = Operator('F')
    G = Operator('G')
