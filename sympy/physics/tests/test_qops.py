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
from sympy.core.add import Add

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
    a*E 
    assert isinstance(a*E, QMul)
    E*F
    assert isinstance(E*F, QMul)
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
        assert issubclass(qmul.evaluates, evaluates)
        assert qmul.hilbert_space == HilbertSpace()

    def should_except(item1, item2):
        try:
            item1*item2
            assert False
        except QuantumError:
            assert True

    a = Symbol('a')
    b = Symbol('b')
    ket = Ket('A')
    ket1 = Ket('J')
    bra = Bra('B')
    bra2 = Bra('l')
    op = Operator('C')
    inner = InnerProduct(bra, ket)
    qmket = op*ket
    qmbra = bra*op
    qmop = op*op
    
    assert op*op == qmket._new_rawargs(Operator, HilbertSpace(), op, op)
    assert op*ket == qmket._new_rawargs(Ket, HilbertSpace(), op, ket)
    should_except(op, bra)
    assert op*inner == qmket._new_rawargs(Operator, HilbertSpace(), op, inner)
    
    should_except(ket, op)
    assert ket*bra == qmket._new_rawargs(Operator, HilbertSpace(), ket, bra)
    should_except(ket, ket)
    assert ket*inner == qmket._new_rawargs(Ket, HilbertSpace(), ket, inner)
    
    should_except(bra, bra)
    assert bra*op == qmket._new_rawargs(Bra, HilbertSpace(), bra, op)
    assert bra*ket == qmket._new_rawargs(InnerProduct, HilbertSpace(), bra, ket)
    assert bra*inner == qmket._new_rawargs(Bra, HilbertSpace(), bra, inner)

    assert inner*op == qmket._new_rawargs(Operator, HilbertSpace(), inner, op)
    assert inner*bra == qmket._new_rawargs(Bra, HilbertSpace(), inner, bra)
    assert inner*ket == qmket._new_rawargs(Ket, HilbertSpace(), inner, ket)
    assert inner*inner == qmket._new_rawargs(InnerProduct, HilbertSpace(), inner, inner)
     
    should_except(qmket, ket)
    helper(qmket*bra, Operator) 
    should_except(qmket, op)
    helper(qmket*inner, Ket)

    helper(qmbra*ket, InnerProduct)
    should_except(qmbra, bra)
    helper(qmbra*op, Bra)
    helper(qmbra*inner, Bra)

    helper(qmop*op, Operator)
    helper(qmop*ket, Ket)
    should_except(qmop, bra)
    helper(qmop*inner, Operator)
    
    assert expand((ket+ket1)*(bra+bra2)) == (ket*bra + ket*bra2 + ket1*bra + ket1*bra2)
    assert Dagger(op*ket) == Dagger(ket)*Dagger(op) == Bra('A')*Dagger(op)

def test_QAdd_mixed():
    def helper(qmul, evaluates):
        assert isinstance(qmul, QAdd)
        assert issubclass(qmul.evaluates, evaluates)
        assert qmul.hilbert_space == HilbertSpace()

    def should_except(item1, item2):
        try:
            item1+item2
            assert False
        except QuantumError:
            assert True

    a = Symbol('a')
    b = Symbol('b')
    ket = Ket('A')
    ket1 = Ket('J')
    bra = Bra('B')
    bra1 = Bra('l')
    op = Operator('C')
    op1 = Operator('d')
    inner = InnerProduct(bra, ket)
    qmop = op + op1
    qmbra = bra + bra1
    qmket = ket + ket1
    possibilities1 = [a, qmket, qmbra, qmop]
    possibilities2 = [inner, op*ket, bra*op, ket*bra]
    assert ket + ket1 == qmop._new_rawargs(Ket, HilbertSpace(), ket, ket1)
    assert bra + bra1 == qmop._new_rawargs(Bra, HilbertSpace(), bra, bra1)
    assert op + op1 == qmop._new_rawargs(Operator, HilbertSpace(), op1, op)
    assert inner + a == Add(a, inner)
    helper(qmket + op*ket, Ket)
    helper(qmbra + possibilities2[2], Bra)
    helper(possibilities1[3] + possibilities2[3], Operator)

    for i in range(4):
        for j in range(4):
            if i != j:
                should_except(possibilities1[i], possibilities2[j])
                should_except(possibilities2[i], possibilities2[j])
    
def test_QPow_mixed():
    def helper(qmul, evaluates, arg1, arg2):
        assert isinstance(qmul, QPow)
        assert issubclass(qmul.evaluates, evaluates)
        assert qmul.base == arg1
        assert qmul.exp == arg2
        assert qmul.hilbert_space == HilbertSpace()

    def should_except(item1, item2):
        try:
            item1**item2
            assert False
        except QuantumError:
            assert True

    a = Symbol('a')
    b = Symbol('b')
    ket = Ket('A')
    ket1 = Ket('J')
    bra = Bra('B')
    bra1 = Bra('l')
    op = Operator('C')
    op1 = Operator('d')
    inner = InnerProduct(bra, ket)
    qmop = op + op1
    qmbra = bra + bra1
    qmket = ket + ket1
    should_except(qmbra, bra)
    should_except(a, bra)
    should_except(bra, op)
    should_except(op, qmbra)
    should_except(qmbra, qmket)
    should_except(qmket, qmbra)
    should_except(a, ket)
    should_except(ket, a)
    should_except(ket, ket1)
    should_except(op, op1)
    should_except(op, ket)
    should_except(ket, op)
    helper(a**op, Operator, a, op)
    helper(op**inner, Operator, op, inner)
    helper(inner**op, Operator, inner, op)
    helper(op**a, Operator, op, a)
       
      
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
