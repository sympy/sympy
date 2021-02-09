from sympy import symbols
from sympy.tensor.array.expressions.array_expressions import ZeroArray
from sympy.testing.pytest import raises


def test_zero_array():
    assert ZeroArray() == 0
    assert ZeroArray().is_Integer

    za = ZeroArray(3, 2, 4)
    assert za.shape == (3, 2, 4)
    za_e = za.as_explicit()
    assert za_e.shape == (3, 2, 4)

    m, n, k = symbols("m n k")
    za = ZeroArray(m, n, k, 2)
    assert za.shape == (m, n, k, 2)
    raises(ValueError, lambda: za.as_explicit())
