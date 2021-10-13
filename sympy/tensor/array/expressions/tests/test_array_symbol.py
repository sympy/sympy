import pytest

from sympy import symbols
from sympy.tensor.array.expressions.array_expressions import ArrayElement, ArraySymbol

k, m, n = symbols("k m n")

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
            )
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
