from sympy.physics.quantum import (
    Operator,
    Dagger
)

from sympy import Matrix, I

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
