from sympy.physics.quantum import (
    Operator,
    Dagger,
    State,
    Ket,
    Bra,
    TimeDepKet,
    TimeDepBra,
    BasisSet
    InnerProduct,
    OuterProduct
)

from sympy import Matrix, I, E, Symbol, expand
from sympy.physics.quantumbasic import QuantumError, QuantumBasic

def test_Dagger():
    i = I
    add = 42+I
    mul = 42*I
    num = 42
    power1 = E**I
    power2 = I**42
    mat = Matrix(((i, add), (mul, num), (power1, power2)))
    A = Bra('A')
    B = Ket('B')
    assert Dagger(i) == -I
    assert Dagger(add) == 42-I
    assert Dagger(mul) == -42*I
    assert Dagger(num) == 42
    assert Dagger(power1) == E**(-I)
    assert Dagger(power2) == -1
    assert Dagger(mat) == Matrix(((-I, -42*I, E**(-I)), (42-I, 42, -1))) == Matrix(((Dagger(i), Dagger(mul), Dagger(power1)), (Dagger(add), Dagger(num), Dagger(power2))))
    assert Dagger(A) == Ket('A')
    assert Dagger(B) == Bra('B')
    assert isinstance(Dagger(A), Ket)
    assert isinstance(Dagger(B), Bra)

def test_Ket():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    assert isinstance(A, Ket)
    assert isinstance(B, Ket)
    assert A.name == Symbol('A')
    assert B.name == Symbol('B')
    assert a*A
    assert a*A == A*a
    assert A+B
    assert A+B == B+A
#    assert expand(a*(A+B)) == a*A+a*B
    assert Dagger(a*A+b*B) == Dagger(A)*Dagger(a)+Dagger(B)*Dagger(b) == Bra('A')*Dagger(a)+Bra('B')*Dagger(b)

def test_Bra():
    a = Symbol('a')
    A = Bra('A')
    B = Bra('B')
    assert isinstance(A, Bra)
    assert isinstance(B, Bra)
    assert A.name == Symbol('A')
    assert B.name == Symbol('B')
    assert a*A
    assert a*A == A*a
    assert A+B
    assert A+B == B+A
#    assert expand(a*(A+B)) == a*A+a*B

def test_InnerProduct():
    A = Bra('A')
    B = Ket('B')
    C = InnerProduct(A, B)
    assert isinstance(C, InnerProduct)
    assert C.bra == A == Bra('A')
    assert C.ket == B == Ket('B')
    assert C != B*A
    assert C == InnerProduct(A, B) == InnerProduct(Bra('A'), Ket('B'))
    assert Dagger(C) == InnerProduct(Dagger(B), Dagger(A)) == InnerProduct(Bra('B'), Ket('A'))
    assert isinstance(Dagger(C), InnerProduct)
    assert C.subs(Bra('A'), Bra('S')) == InnerProduct(Bra('S'), Ket('B')) == InnerProduct(Bra('S'), B)
    assert C.subs(Ket('B'), Ket('S')) == InnerProduct(Bra('A'), Ket('S')) == InnerProduct(A, Ket('S'))

def test_OuterProduct():
    A = Ket('A')
    B = Bra('B')
    C = OuterProduct(A, B)
    assert isinstance(C, OuterProduct)
    assert C.bra == B == Bra('B')
    assert C.ket == A == Ket('A')
    assert C != B*A
    assert C == OuterProduct(A, B) == OuterProduct(Ket('A'), Bra('B'))
    assert Dagger(C) == OuterProduct(Dagger(B), Dagger(A)) == OuterProduct(Ket('B'), Bra('A'))
    assert isinstance(Dagger(C), OuterProduct)
    assert C.subs(Bra('B'), Bra('S')) == OuterProduct(Ket('A'), Bra('S')) == OuterProduct(A, Bra('S'))
    assert C.subs(Ket('A'), Ket('S')) == OuterProduct(Ket('S'), Bra('B')) == OuterProduct(Ket('S'), B)

def test_Operator():
    a = Symbol('a')
    A = Operator('A')
    B = Operator('B')
    C = Operator('C')
    assert isinstance(A, Operator)
    assert isinstance(B, Operator)
    assert a*A
    assert A+B == B+A
    assert (A+B)+C == A+(B+C)
#    assert expand(a*(A+B)) == a*A+a*B
    assert Dagger(a*A) == Dagger(A)*Dagger(a)
    assert A*B != B*A
    assert A*(B*C) == (A*B)*C

def test_mixed():
    A = Ket('A')
    B = Ket('B')
    C = Operator('C')
    D = Operator('D')
    assert C*B
#    assert expand((C+D)*(A+B)) == (C*A + C*B + D*A + D*B)
#    assert Dagger(C*A) == Dagger(A)*Dagger(C) == Bra('A')*Dagger(C)
