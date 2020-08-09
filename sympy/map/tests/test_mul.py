from sympy import (
    AdditionOperator, MultiplicationOperator, Set, S, symbols,
    NumericMultiplicationOperator, VectorAdditionOperator,
    ScalarMultiplicationOperator, AbelianGroup,
)

def test_MultiplicationOperator():
    A = Set('A')
    a, b, c, e1, e2 = [A.element(n) for n in ('a', 'b', 'c', 'e1', 'e2')]
    add = AdditionOperator(A**2, A, e1)
    mul = MultiplicationOperator(A**2, A, e2)
    pow = mul.exponent_operator()
    div = mul.right_division_operator()

    # evaluating commutative multiplication sorts argument
    assert mul(c, b, a,  add_op=add, evaluate=True).arguments == (a, b, c)

    # Nested structure is flattened
    assert mul(a, mul(b, c, add_op=add),  add_op=add, evaluate=True).arguments == (a, b, c)
    # Identity is removed
    assert mul(a, e2, b, add_op=add, evaluate=True).arguments == (a, b)

    # Generalized power
    assert mul(a, a, add_op=add, evaluate=True) == pow(a, 2)

    # as_coeff_Mul
    assert mul(a, b, add_op=add).as_coeff_Mul() == (e2, mul(a, b, add_op=add))

    # distribution
    assert mul.distribute([a, add(b, c, mul_op=mul)], add_op=add) \
        == add(mul(a, b, add_op=add), mul(a, c, add_op=add), mul_op=mul)

    # inverse
    assert div(b, e2, add_op=add, evaluate=True) == b
    assert div(e2, b, add_op=add, evaluate=True) == pow(b, -1)
    assert div(a, b, add_op=add, evaluate=True) == mul(a, pow(b, -1), add_op=add, evaluate=True)
    assert div(a, a, add_op=add, evaluate=True) == e2

def test_NumericMultiplicationOperator():
    x, y = symbols('x y', real=True)
    mul = NumericMultiplicationOperator(S.Reals**2, S.Reals)
    pow = mul.exponent_operator()
    div = mul.right_division_operator()

    # normal multiplication
    assert mul(1, 2, 3, 4, evaluate=True) == 24
    # numbers come first
    assert mul(x, 2, y, 4, evaluate=True).arguments == (8, x, y)
    # Identity is removed
    assert mul(x, 1, y, evaluate=True).arguments == (x, y)

    # as_coeff_mul
    assert mul(x, y).as_coeff_Mul() == (1, mul(x,y))
    assert mul(x, 3).as_coeff_Mul() == (1, mul(x, 3))
    assert mul(x, 3, evaluate=True).as_coeff_Mul() == (3, x)

    # inverse
    assert div(1, 3, evaluate=True) == pow(3, -1)
    assert div(3, 1, evaluate=True) == mul(3, 1, evaluate=True) == 3
    assert pow(1, -1, evaluate=True) == pow(1, x, evaluate=True) == 1
    assert div(1, 1, evaluate=True) == 1

def test_ScalarMultiplicationOperator():
    V = Set('V')
    v, V_e = [V.element(i) for i in ('v', 'e')]
    vv_add = VectorAdditionOperator(V**2, V, V_e)
    V_group = AbelianGroup('V', (V,), (vv_add,))
    mul = ScalarMultiplicationOperator(S.RealsField*V_group, V_group)
    assert mul(2, mul(3, v), evaluate=True) == mul(6, v)
