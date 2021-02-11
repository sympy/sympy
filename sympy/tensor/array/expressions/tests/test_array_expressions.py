from sympy import symbols
from sympy.tensor.array.expressions.array_expressions import ZeroArray, OneArray
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


def test_one_array():
    assert OneArray() == 1
    assert OneArray().is_Integer

    oa = OneArray(3, 2, 4)
    assert oa.shape == (3, 2, 4)
    oa_e = oa.as_explicit()
    assert oa_e.shape == (3, 2, 4)

    m, n, k = symbols("m n k")
    oa = OneArray(m, n, k, 2)
    assert oa.shape == (m, n, k, 2)
    raises(ValueError, lambda: oa.as_explicit())
