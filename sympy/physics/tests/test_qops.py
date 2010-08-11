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
from sympy.physics.qmul import QMul
from sympy.physics.qadd import QAdd
from sympy.physics.qpow import QPow
from sympy import Symbol, expand
from sympy.utilities.pytest import XFAIL

def test_QMul_Ket():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    assert a*A
    assert isinstance(a*A, QMul)
    assert (a*A).evaluates == A.__class__ == Ket
    assert a*A == A*a
    assert expand(a*(A+B)) == a*A+a*B

@XFAIL
def test_QMul_Ket_fail():
    A = Ket('a')
    B = Ket('b')
    assert A*B

def test_QMul_Bra():
    a = Symbol('a')
    b = Symbol('b')
    C = Bra('C')
    D = Bra('D')
    assert a*C
    assert isinstance(a*C, QMul)
    assert (a*C).evaluates == C.__class__ == Bra
    assert a*C == C*a
    assert expand(a*(C+D)) == a*C+a*D

@XFAIL
def test_QMul_Bra_fail():
    A = Bra('a')
    B = Bra('b')
    assert A*B

def test_QMul_Operator():
    a = Symbol('a')
    b = Symbol('b')
    E = Operator('E')
    F = Operator('F')
    G = Operator('G')
    assert a*E
    assert (a*E, QMul)
    assert E*F
    assert (E*F, QMul)
    assert (E*F).evaluates == E.__class__ == Operator
    assert a*E == E*a
    assert expand(a*(E+F)) == a*E+a*F
    assert E*F != F*E
    assert E*(F*G) == (E*F)*G
    assert Dagger(a*E) == Dagger(E)*Dagger(a)
    assert Dagger(a*E).evaluates == E.__class__ == Operator

def test_QMul_InnerProduct():
    inner1 = InnerProduct(Bra('a'), Ket('b'))
    inner2 = InnerProduct(Bra('c'), Ket('d'))
    prod = inner1*inner2
    assert prod*inner1
    assert isinstance(prod, QMul)
    assert prod.evaluates == inner1.__class__ == InnerProduct

def test_QMul_OuterProduct():
    outer1 = OuterProduct(Ket('a'), Bra('b'))
    outer2 = OuterProduct(Ket('c'), Bra('d'))
    prod = outer1*outer2
    assert prod*outer1
    assert isinstance(prod, QMul)
    assert prod.evaluates == outer1.__class__ == OuterProduct

def test_QMul_mixed():
    def helper(qmul, evaluates):
        assert isinstance(qmul, QMul)
        assert qmul.evaluates == evaluates
        assert qmul.hilbert_space = HilbertSpace()
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    C = Bra('C')
    E = Operator('E')
    F = Operator('F')
    G = Operator('G')
    H = InnerProduct(C, B)
    I = InnerProduct(C, A)
    qmket = E*B
    qmbra = C*F
    qmop = F*G
    assert E*B == qmket._new_rawargs(Ket, HilbertSpace(), 1, E, B)
    assert C*F == qmket._new_rawargs(Bra, HilbertSpace(), 1, C, F)
    assert C*F == qmket._new_rawargs(Bra, HilbertSpace(), 1, C, F)
    assert A*H == qmket._new_rawargs(Ket, HilbertSpace(), 1, A, H)
    assert C*H == qmket._new_rawargs(Bra, HilbertSpace(), 1, C, H)
    assert E*H == qmket._new_rawargs(Ket, HilbertSpace(), 1, E, H)
    assert I*H == qmket._new_rawargs(InnerProduct, HilbertSpace(), 1, I, H)
    assert H*A == qmket._new_rawargs(Ket, HilbertSpace(), 1, H, A)
    assert H*C == qmket._new_rawargs(Bra, HilbertSpace(), 1, H, C)
    assert H*E == qmket._new_rawargs(Operator, HilbertSpace(), 1, H, E)
    assert H*I == qmket._new_rawargs(InnerProduct, HilbertSpace(), 1, H, I)
    helper(E*qmket, Ket)
    assert qmbra*F == qmket._new_rawargs(Bra, HilbertSpace(), 1, qmbra, F)
    assert E*qmop == qmket._new_rawargs(Operator, HilbertSpace(), 1, E, qmop)
    assert expand((E+F)*(A+B)) == (E*A + E*B + F*A + F*B)
    assert Dagger(E*A) == Dagger(A)*Dagger(E) == Bra('A')*Dagger(E)

@XFAIL
def test_QMul_mixed_fail():
    A = Ket('a')
    B = Bra('b')
    C = Operator
    assert A*C
    assert C*B

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
    assert E**a
    assert a**F

"""
    assert A*C*E*H
    assert isinstance(A*C*E*H, QMul)
    assert A*C*H*E
    assert A*H*E*C
    assert A*H*C*E
    assert A*E*C*H
    assert A*E*H*C
    assert H*A*E*C
    assert H*A*C*E
    assert H*C*E*A
    assert H*C*A*E
    assert H*E*C*A
    assert H*E*A*C
    assert C*A*H*E
    assert C*A*E*H
    assert C*H*A*E
    assert C*H*E*A
    assert C*E*A*H
    assert C*E*H*A
    assert E*A*C*H
    assert E*A*H*C
    assert E*C*A*H
    assert E*C*H*A
    assert E*H*C*A
    assert E*H*A*C
"""
