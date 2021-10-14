from typing import Type

import pytest

from sympy import Tuple, symbols
from sympy.tensor.array.expressions.array_expressions import (
    ArrayElement,
    ArraySlice,
    ArraySymbol,
)

k, m, n = symbols("k m n")


class TestArrayElement:
    @pytest.mark.parametrize(
        ["shape", "indices", "expected"],
        [
            # shapeless
            ([], (0,), (0,)),
            ([], (1, 2), (1, 2)),
            # specific shape
            ([3, 3], (1, 2), (1, 2)),
            ([3, 3, 3], (1, 2), (1, 2)),
            # negative indices
            ([], (-1,), (-1,)),
            ([3], (-1,), (2,)),
            ([3, 3], (-1, -2), (2, 1)),
        ],
    )
    def test_construct_from_array_symbol(self, shape, indices, expected):
        A = ArraySymbol("A", *shape)
        element = ArrayElement(A, indices=indices)
        assert element.parent is A
        assert element.indices == expected

    @pytest.mark.parametrize(
        ["shape", "indices", "exception", "match"],
        [
            ([2], (3,), ValueError, "shape is out of bounds"),
            ([3], (0, 0), IndexError, r"Too many indices for ArrayElement"),
        ],
    )
    def test_construction_errors_from_array_symbol(
        self, shape, indices, exception: Type[Exception], match
    ):
        A = ArraySymbol("A", *shape)
        with pytest.raises(exception, match=match):
            ArrayElement(A, indices=indices)


class TestArraySlice:
    @pytest.mark.parametrize(
        ["shape", "indices", "expected"],
        [
            # shapeless
            ([], (0,), (0,)),
            ([], (slice(None),), (Tuple(0, None, None),)),
            ([], (slice(3, None),), (Tuple(3, None, None),)),
            ([], (slice(3, None, 2),), (Tuple(3, None, 2),)),
            ([], (slice(None, 5),), (Tuple(0, 5, None),)),
            ([], (1, slice(None, 5)), (1, Tuple(0, 5, None))),
            # # specific shape (normalized)
            ([3, 3], (1, slice(None, 2)), (1, Tuple(0, 2, 1))),
            ([3, 3], (1, slice(None, 5)), (1, Tuple(0, 5, 1))),  # overflow
            # # negative indices
            ([3, 3], (-1, slice(None, -2)), (-1, Tuple(0, 1, 1))),
        ],
    )
    def test_construct_from_array_symbol(self, shape, indices, expected):
        A = ArraySymbol("A", *shape)
        array_slice = ArraySlice(A, slices=indices)
        assert array_slice.parent is A
        assert array_slice.slices == expected


class TestArraySymbol:
    @pytest.mark.parametrize(
        "shape",
        [
            (),
            (3, 2, 4),
            (k, m, n),
        ],
    )
    def test_constructor(self, shape: tuple):
        A = ArraySymbol("A", *shape)
        assert A.name == "A"
        assert A.shape == shape

    def test_equality(self):
        assert ArraySymbol("A") == ArraySymbol("A")
        assert ArraySymbol("A") != ArraySymbol("B")
        assert ArraySymbol("A", 2, 3) != ArraySymbol("A")
        assert ArraySymbol("A", m) != ArraySymbol("A", n)
        assert ArraySymbol("A", n) == ArraySymbol("A", n)

    def test_getitem(self):
        A = ArraySymbol("A")
        assert A[0] == ArrayElement(A, indices=(0,))
        assert A[9, 7] == ArrayElement(A, indices=(9, 7))
        A = ArraySymbol("A", n)
        assert A[2] == ArrayElement(A, indices=(2,))
        A = ArraySymbol("A", 3, 2, 4)
        assert A[2, 1] == ArrayElement(A, indices=(2, 1))
        assert A[2, 1, 3] == ArrayElement(A, indices=(2, 1, 3))
        assert A[-2, -1, 3] == ArrayElement(A, indices=(1, 1, 3))
        with pytest.raises(ValueError, match="shape is out of bounds"):
            A[4, 1, 3]
        with pytest.raises(
            IndexError,
            match=(
                f"Too many indices for {ArrayElement.__name__}: parent "
                f"{ArraySymbol.__name__} is 3-dimensional, but 4 indices were given"
            ),
        ):
            A[0, 1, 2, 5]

    def test_getitem_slice(self):
        A_slice = ArraySymbol("A")[1:3, :5]
        assert A_slice.shape == (2, 5)
        A_slice = ArraySymbol("A")[1:n, m, 3:7:2]
        assert A_slice.shape == (n - 1, 1, 2)
        A_slice = ArraySymbol("A")[3:]
        assert A_slice.shape == (None,)

    def test_getitem_slice_overflow(self):
        assert ArraySymbol("A", 3)[-5:].shape == (3,)
        assert ArraySymbol("A", 3)[:5].shape == (3,)
        assert ArraySymbol("A", 3)[:5:2].shape == (2,)
        assert ArraySymbol("A", n)[-5:].shape == (5,)
        assert ArraySymbol("A", n)[:5].shape == (5,)
        assert ArraySymbol("A", n)[:5:2].shape == (2,)
        assert ArraySymbol("A", n, n)[:5:2].shape == (2, n)
