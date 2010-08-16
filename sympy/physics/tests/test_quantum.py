from sympy.physics.quantum import (
    Operator,
    Dagger,
    State,
    Ket,
    Bra,
    TimeDepState,
    TimeDepKet,
    TimeDepBra,
    BasisSet,
    InnerProduct,
    OuterProduct,
    KroneckerDelta,
    Commutator,
    split_commutative_parts
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
from sympy import Matrix, I, E, Symbol, expand, Number, symbols
from sympy.physics.quantumbasic import QuantumError, QuantumBasic

def test_Dagger():
    i = I
    add = 42+I
    mul = 42*I
    num = 42
    power1 = E**I
    power2 = I**42
    mat = Matrix(((i, add), (mul, num), (power1, power2)))
    assert Dagger(i) == -I
    assert Dagger(add) == 42-I
    assert Dagger(mul) == -42*I
    assert Dagger(num) == 42
    assert Dagger(power1) == E**(-I)
    assert Dagger(power2) == -1
    assert Dagger(mat) == Matrix(((-I, -42*I, E**(-I)), (42-I, 42, -1))) ==\
    Matrix(((Dagger(i), Dagger(mul), Dagger(power1)), (Dagger(add),\
    Dagger(num), Dagger(power2))))

def test_Ket():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    assert isinstance(A, Ket)
    assert A.name == Symbol('A')
    assert A.is_symbolic == True
    assert isinstance(A.dual, Bra)
    assert A.dual == Bra('A')
    assert A.hilbert_space == HilbertSpace()
    assert A.acts_like == A.__class__ == Ket
    assert Dagger(B) == Bra('B')
    assert isinstance(Dagger(B), Bra)
    assert Dagger(B).hilbert_space == HilbertSpace()
    assert Dagger(a*A+b*B) == Dagger(A)*Dagger(a)+Dagger(B)*Dagger(b) ==\
    Bra('A')*Dagger(a)+Bra('B')*Dagger(b)


def test_TimeDepKet():
    A = TimeDepKet('A', 't')
    assert isinstance(A, TimeDepKet) == isinstance(A, TimeDepState)
    assert A.name == Symbol('A')
    assert A.hilbert_space == HilbertSpace()
    assert A.is_symbolic == True
    assert A.acts_like == A.__class__ == TimeDepKet
    assert isinstance(A.dual, TimeDepBra)
    assert A.dual == TimeDepBra('A', 't')

def test_Bra():
    a = Symbol('a')
    b = Symbol('b')
    A = Bra('A')
    B = Bra('B')
    assert isinstance(A, Bra)
    assert A.name == Symbol('A')
    assert A.is_symbolic == True
    assert isinstance(A.dual, Ket)
    assert A.dual == Ket('A')
    assert A.hilbert_space == HilbertSpace()
    assert A.acts_like == A.__class__ == Bra
    assert Dagger(A) == Ket('A')
    assert isinstance(Dagger(A), Ket)
    assert Dagger(A).hilbert_space == HilbertSpace()
    assert Dagger(a*A+b*B) == Dagger(A)*Dagger(a)+Dagger(B)*Dagger(b) ==\
    Ket('A')*Dagger(a)+Ket('B')*Dagger(b)


def test_TimeDepBra():
    A = TimeDepBra('A', 't')
    assert isinstance(A, TimeDepBra) == isinstance(A, TimeDepState)
    assert A.name == Symbol('A')
    assert A.hilbert_space == HilbertSpace()
    assert A.is_symbolic == True
    assert A.acts_like == A.__class__ == TimeDepBra
    assert isinstance(A.dual, TimeDepKet)
    assert A.dual == TimeDepKet('A', 't')

def test_InnerProduct():
    A = Bra('A')
    B = Ket('B')
    C = InnerProduct(A, B)
    assert isinstance(C, InnerProduct)
    assert C.bra == A == Bra('A')
    assert C.ket == B == Ket('B')
    assert C != B*A
    assert C == InnerProduct(Bra('A'), Ket('B'))
    assert C.acts_like == C.__class__ == InnerProduct
    assert C.hilbert_space == HilbertSpace()
    assert Dagger(C) == InnerProduct(Dagger(B), Dagger(A)) ==\
    InnerProduct(Bra('B'), Ket('A'))
    assert isinstance(Dagger(C), InnerProduct)
    assert C.subs(Bra('A'), Bra('S')) == InnerProduct(Bra('S'), Ket('B')) ==\
    InnerProduct(Bra('S'), B)
    assert C.subs(Ket('B'), Ket('S')) == InnerProduct(Bra('A'), Ket('S')) ==\
    InnerProduct(A, Ket('S'))

def test_OuterProduct():
    A = Ket('A')
    B = Bra('B')
    C = OuterProduct(A, B)
    assert isinstance(C, OuterProduct) == isinstance(C, Operator)
    assert C.bra == B == Bra('B')
    assert C.ket == A == Ket('A')
    assert C != B*A
    assert C == OuterProduct(Ket('A'), Bra('B'))
    assert C.acts_like == C.__class__ == OuterProduct
    assert C.hilbert_space == HilbertSpace()
    assert Dagger(C) == OuterProduct(Dagger(B), Dagger(A)) ==\
    OuterProduct(Ket('B'), Bra('A'))
    assert isinstance(Dagger(C), OuterProduct)
    assert C.subs(Bra('B'), Bra('S')) == OuterProduct(Ket('A'), Bra('S')) ==\
    OuterProduct(A, Bra('S'))
    assert C.subs(Ket('A'), Ket('S')) == OuterProduct(Ket('S'), Bra('B')) ==\
    OuterProduct(Ket('S'), B)

def test_Operator():
    a = Symbol('a')
    A = Operator('A')
    B = Operator('B')
    C = Operator('C')
    D = Bra('D')
    assert isinstance(A, Operator)
    assert A.name == Symbol('A')
    assert A.acts_like == A.__class__ == Operator
    assert A.hilbert_space == HilbertSpace()
    assert isinstance(Dagger(C), Dagger)
    assert Dagger(C).acts_like == C.__class__ == Operator
    assert Dagger(C).subs(C, D) == Dagger(D) == Ket('D')
    assert isinstance(Dagger(Dagger(C)), Operator)

def test_KroneckerDelta():
    assert KroneckerDelta(1,2) == 0
    assert KroneckerDelta(1,1) == 1
    assert Dagger(KroneckerDelta(1,2)) == 0
    assert Dagger(KroneckerDelta(1,1)) == 1

def test_Commutator():
    A = Symbol('A', **{'commutative':False})
    B = Symbol('B', **{'commutative':False})
    C = Symbol('C', **{'commutative':False})
    x,y = symbols('xy')
    bra = Bra('b')
    ket = Ket('a')
    E = Operator('E')
    F = Operator('F')
    assert Commutator(0, 2) == 0
    assert Commutator(A+B, C) == Commutator(A, C) + Commutator(B, C)
    assert Commutator(x*A, y*B) == x*y*Commutator(A, B)
