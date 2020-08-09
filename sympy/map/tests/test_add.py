from sympy import (
    AdditionOperator, MultiplicationOperator, Set, S, symbols,
    NumericAdditionOperator, NumericMultiplicationOperator
)

def test_AdditionOperator():
    A = Set('A')
    a, b, c, e1, e2 = [A.element(n) for n in ('a', 'b', 'c', 'e1', 'e2')]
    add = AdditionOperator(A**2, A, e1)
    repeat_add = add.exponent_operator()
    add_inv = add.inverse_operator()
    sub = add.subtraction_operator()
    mul = MultiplicationOperator(A**2, A, e2)

    # Nested structure is flattened
    assert add(a, add(b, c, mul_op=mul),  mul_op=mul, evaluate=True).arguments == (a, b, c)
    # Identity is removed
    assert add(a, e1, b, mul_op=mul, evaluate=True).arguments == (a, b)

    # Un-distribution
    assert add(a, a, mul_op=mul, evaluate=True) == mul(add(e2, e2, mul_op=mul), a, add_op=add)

    # as_coeff_Add
    assert add(a, b, mul_op=mul).as_coeff_Add() == (e1, add(a, b, mul_op=mul))

    # inverse
    assert sub(a, e1, mul_op=mul, evaluate=True) == a
    assert sub(e1, a, mul_op=mul, evaluate=True) == add_inv(a)
    assert sub(a, b, mul_op=mul) == add(a, add_inv(b), mul_op=mul)
    assert sub(a, add_inv(b), mul_op=mul, evaluate=True) == add(a, b, mul_op=mul)
    assert sub(a, a, mul_op=mul, evaluate=True) == e1
    assert sub(repeat_add(a, 2), a, mul_op=mul, evaluate=True) == a

def test_NumericAdditionOperator():
    x, y = symbols('x y', real=True)
    add = NumericAdditionOperator(S.Reals**2, S.Reals)
    add_inv = add.inverse_operator()
    sub = add.subtraction_operator()
    mul = NumericMultiplicationOperator(add.domain, add.codomain)

    # normal addition
    assert add(1, 2, 3, 4, evaluate=True) == 10
    # numbers come first
    assert add(x, 2, y, 4, evaluate=True).arguments == (6, x, y)
    # Identity is removed
    assert add(x, 0, y, evaluate=True).arguments == (x, y)

    # undistribution
    assert add(x, x, evaluate=True) == mul(2, x)

    # as_coeff_Add
    assert add(x, y).as_coeff_Add() == (0, add(x,y))
    assert add(x, 3).as_coeff_Add() == (0, add(x, 3))
    assert add(x, 3, evaluate=True).as_coeff_Add() == (3, x)

    # inverse
    assert add_inv(3, evaluate=True) == -3
    assert sub(1, 2) == add(1, -2)
    assert sub(1, 2, evaluate=True) == -1
