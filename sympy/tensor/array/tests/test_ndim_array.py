from __future__ import annotations
from sympy.testing.pytest import raises
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.matrices.dense import Matrix
from sympy.simplify import simplify
from sympy.tensor.array import Array
from sympy.tensor.array.dense_ndim_array import (
    ImmutableDenseNDimArray, MutableDenseNDimArray)
from sympy.tensor.array.sparse_ndim_array import (
    ImmutableSparseNDimArray, MutableSparseNDimArray)

from sympy.abc import x, y

mutable_array_types = [
    MutableDenseNDimArray,
    MutableSparseNDimArray
]

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


def test_issue_20222():
    A = Array([[1, 2], [3, 4]])
    B = Matrix([[1,2],[3,4]])
    raises(TypeError, lambda: A - B)


def test_issue_17851():
    for array_type in array_types:
        A = array_type([])
        assert isinstance(A, array_type)
        assert A.shape == (0,)
        assert list(A) == []


def test_issue_and_18715():
    for array_type in mutable_array_types:
        A = array_type([0, 1, 2])
        A[0] += 5
        assert A[0] == 5


def test_applyfunc_scalar_and_iterator_input():
    # applyfunc on a 0-dim (scalar) array used to crash: it builds the new
    # array from a ``map`` object, and the shape==() length check called
    # ``len()`` on that one-shot iterator.  This also covers passing a bare
    # iterator (map/generator) together with an explicit shape.
    for array_type in array_types:
        scalar = array_type(5)
        assert scalar.shape == ()
        assert array_type(5).applyfunc(lambda e: e + 1)[()] == 6
        assert scalar.applyfunc(lambda e: e**2)[()] == 25
        # a one-shot iterator with an explicit shape must be accepted:
        assert array_type(map(lambda e: e, [7]), ()) [()] == 7
        assert list(array_type((i for i in range(3)), (3,))) == [0, 1, 2]
