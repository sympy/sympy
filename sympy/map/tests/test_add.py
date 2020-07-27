from sympy import S, symbols
from sympy.map import scalar_add, scalar_mul
x, y, z = symbols('x y z')

def test_AdditionOperator():
    # If mul_op is not given, repetitive elements are not converted.
    assert scalar_add(x, x, evaluate=True).arguments == (x, x)
    assert scalar_add(x, x, mul_op=scalar_mul, evaluate=True) == scalar_mul(2, x)

    # Nested addition is always flattened
    assert scalar_add(x, scalar_add(y, z), evaluate=True).arguments == (x, y, z)
    # Identity is always removed
    assert scalar_add(x, 0, y, evaluate=True).arguments == (x, y)
