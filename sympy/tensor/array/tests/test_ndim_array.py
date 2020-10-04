from sympy.testing.pytest import raises
from sympy import (
    Array, ImmutableDenseNDimArray, ImmutableSparseNDimArray,
    MutableDenseNDimArray, MutableSparseNDimArray, sin, cos,
    simplify
)
from sympy.abc import x, y

array_types = [
    ImmutableDenseNDimArray,
    ImmutableSparseNDimArray,
    MutableDenseNDimArray,
    MutableSparseNDimArray
]


def test_array_negative_indices():
    for ArrayType in array_types:
        test_array = ArrayType([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
        assert test_array[:, -1] == Array([5, 10])
        assert test_array[:, -2] == Array([4, 9])
        assert test_array[:, -3] == Array([3, 8])
        assert test_array[:, -4] == Array([2, 7])
        assert test_array[:, -5] == Array([1, 6])
        assert test_array[:, 0] == Array([1, 6])
        assert test_array[:, 1] == Array([2, 7])
        assert test_array[:, 2] == Array([3, 8])
        assert test_array[:, 3] == Array([4, 9])
        assert test_array[:, 4] == Array([5, 10])

        raises(ValueError, lambda: test_array[:, -6])
        raises(ValueError, lambda: test_array[-3, :])

        assert test_array[-1, -1] == 10


def test_issue_18361():
    A = Array([sin(2 * x) - 2 * sin(x) * cos(x)])
    B = Array([sin(x)**2 + cos(x)**2, 0])
    C = Array([(x + x**2)/(x*sin(y)**2 + x*cos(y)**2), 2*sin(x)*cos(x)])
    assert simplify(A) == Array([0])
    assert simplify(B) == Array([1, 0])
    assert simplify(C) == Array([x + 1, sin(2*x)])
