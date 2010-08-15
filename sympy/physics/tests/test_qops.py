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
    assert (a*A).acts_like == A.__class__ == Ket
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
    assert (a*C).acts_like == C.__class__ == Bra
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
    assert (E*F).acts_like == E.__class__ == Operator
    assert a*E == E*a
    assert expand(a*(E+F)) == a*E+a*F
    assert E*F != F*E
    assert E*(F*G) == (E*F)*G
    assert Dagger(a*E) == Dagger(E)*Dagger(a)
    assert Dagger(a*E).acts_like == E.__class__ == Operator

def test_QMul_InnerProduct():
    inner1 = InnerProduct(Bra('a'), Ket('b'))
    inner2 = InnerProduct(Bra('c'), Ket('d'))
    prod = inner1*inner2
    assert prod*inner1
    assert isinstance(prod, QMul)
    assert prod.acts_like == inner1.__class__ == InnerProduct

def test_QMul_OuterProduct():
    outer1 = OuterProduct(Ket('a'), Bra('b'))
    outer2 = OuterProduct(Ket('c'), Bra('d'))
    prod = outer1*outer2
    assert prod*outer1
    assert isinstance(prod, QMul)
    assert prod.acts_like == outer1.__class__ == OuterProduct

def test_QMul_mixed():
    def helper(qmul, acts_like):
        assert isinstance(qmul, QMul)
        assert issubclass(qmul.acts_like, acts_like)
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

    assert op*op == op**2
    assert op*ket == QMul._new_rawargs(Ket, HilbertSpace(), op, ket)
    should_except(op, bra)
    assert op*inner == QMul._new_rawargs(Operator, HilbertSpace(), op, inner)

    should_except(ket, op)
    assert ket*bra == QMul._new_rawargs(Operator, HilbertSpace(), ket, bra)
    should_except(ket, ket)
    assert ket*inner == QMul._new_rawargs(Ket, HilbertSpace(), ket, inner)

    should_except(bra, bra)
    assert bra*op == QMul._new_rawargs(Bra, HilbertSpace(), bra, op)
    assert bra*ket == QMul._new_rawargs(InnerProduct, HilbertSpace(), bra, ket)
    assert bra*inner == QMul._new_rawargs(Bra, HilbertSpace(), bra, inner)

    assert inner*op == QMul._new_rawargs(Operator, HilbertSpace(), inner, op)
    assert inner*bra == QMul._new_rawargs(Bra, HilbertSpace(), inner, bra)
    assert inner*ket == QMul._new_rawargs(Ket, HilbertSpace(), inner, ket)
    assert inner*inner == inner**2

    should_except(qmket, ket)
    helper(qmket*bra, Operator)
    should_except(qmket, op)
    helper(qmket*inner, Ket)

    helper(qmbra*ket, InnerProduct)
    should_except(qmbra, bra)
    helper(qmbra*op, Bra)
    helper(qmbra*inner, Bra)

    assert qmop*op == op**3
    helper(qmop*ket, Ket)
    should_except(qmop, bra)
    helper(qmop*inner, Operator)

    assert expand((ket+ket1)*(bra+bra2)) == (ket*bra + ket*bra2 + ket1*bra +\
    ket1*bra2)
    assert Dagger(op*ket) == Dagger(ket)*Dagger(op) == Bra('A')*Dagger(op)

def test_QAdd_mixed():
    def helper(qmul, acts_like):
        assert isinstance(qmul, QAdd)
        assert issubclass(qmul.acts_like, acts_like)
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
    assert ket + ket1 == QAdd._new_rawargs(Ket, HilbertSpace(), ket, ket1)
    assert bra + bra1 == QAdd._new_rawargs(Bra, HilbertSpace(), bra, bra1)
    assert op + op1 == QAdd._new_rawargs(Operator, HilbertSpace(), op1, op)
    assert inner + a == Add(a, inner)
    helper(qmket + op*ket, Ket)
    helper(qmbra + possibilities2[2], Bra)
    helper(possibilities1[3] + possibilities2[3], Operator)

    for i in range(4):
        for j in range(4):
            if i != j:
                should_except(possibilities1[i], possibilities2[j])
                should_except(possibilities2[i], possibilities2[j])

def test_QPow_combinations():
    op = Operator('A')
    op1 = Operator('B')
    ket = Ket('C')
    qmop = op*op1
    assert op*op*op == op**3
    assert op*qmop*op1 == op**2*op1**2
    assert (op*op1)*(op1*ket) == op*op1**2*ket
    assert 2**op*2**op1 != 2**(op+op1)
    assert 2**(op*op1)*2**(op*op) != 2**(op*op1+op**2)

def test_QPow_mixed():
    def helper(qmul, acts_like, arg1, arg2):
        assert isinstance(qmul, QPow)
        assert issubclass(qmul.acts_like, acts_like)
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
