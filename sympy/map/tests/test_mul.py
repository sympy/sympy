from sympy import S, symbols
from sympy.map import scalar_add, scalar_mul
x, y, z = symbols('x y z')

def test_MultiplicationOperator():
    # Nested is always flattened
    assert scalar_mul(x, scalar_mul(y, z), evaluate=True).arguments == (x, y, z)
    # Identity is removed
    assert scalar_mul(x, 1, y, evaluate=True).arguments == (x, y)

    # numbers come first
    assert scalar_mul(x, 2, y, 4, evaluate=True).arguments == (8, x, y)

    # can be distributed
    assert scalar_mul.distribute(
        [x, scalar_add(2, y)], add_op= scalar_add, evaluate=True
    ) == scalar_add(
        scalar_mul(2, x), scalar_mul(x, y)
    )

    # as_coeff_Mul
    assert scalar_mul(2, x, y, evaluate=True).as_coeff_Mul(mul_op=scalar_mul) == (2, scalar_mul(x, y))
