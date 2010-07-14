from sympy.physics.quantum import (
    Operator,
    Dagger,
    State,
    Ket,
    Bra,
    InnerProduct,
    OuterProduct
)

from sympy import Matrix, I, Symbol, expand

def test_Dagger():
    i = I
    add = 42+I
    mul = 42*I
    num = 42
    power1 = 42**I
    power2 = I**42
    mat = Matrix(((i, add), (mul, num), (power1, power2)))
    assert Dagger(i) == -I
    assert Dagger(add) == 42-I
    assert Dagger(mul) == -42*I
    assert Dagger(num) == 42
    assert Dagger(power1) == 42**(-I)
    assert Dagger(power2) == -1
    assert Dagger(mat) == Matrix(((-I, -42*I, 42**(-I)), (42-I, 42, -1)))# == Matrix(((Dagger(i), Dagger(mul), Dagger(power1)), (Dagger(add), Dagger(num), Dagger(power2)))

def test_Ket():
    a = Symbol('a')
    A = Ket('A')
    B = Ket('B')
    assert isinstance(A, Ket)
    assert isinstance(B, Ket)
    assert a*A
    assert a*A == A*a
#    assert A*B != B*A
    assert A+B
    assert A+B == B+A
    assert expand(a*(A+B)) == a*A+a*B

def test_Dagger_Ket():
    a = Symbol('a')
    b = Symbol('b')
    A = Ket('A')
    B = Ket('B')
    assert Dagger(a*A+b*B) == Dagger(A)*Dagger(a)+Dagger(B)*Dagger(b) == Bra('A')*Dagger(a)+Bra('B')*Dagger(b)

def test_Bra():
    a = Symbol('a')
    A = Bra('A')
    B = Bra('B')
    assert isinstance(A, Bra)
    assert isinstance(B, Bra)
    assert a*A
    assert a*A == A*a
#    assert A*B != B*A
    assert A+B
    assert A+B == B+A
    assert expand(a*(A+B)) == a*A+a*B

def test_State():
    a = Symbol('a')
    b = Symbol('b')

def test_InnerProduct():
    A = Bra('A')
    B = Ket('B')
    assert isinstance(A*B, InnerProduct)
    assert A*B != B*A
    assert A*B == InnerProduct(A, B) == InnerProduct(Bra('A'), Ket('B'))

def test_OuterProduct():
    A = Ket('A')
    B = Bra('B')
    assert isinstance(A*B, OuterProduct)
    assert A*B != B*A
    assert A*B == OuterProduct(A, B) == OuterProduct(Ket('A'), Bra('B'))
