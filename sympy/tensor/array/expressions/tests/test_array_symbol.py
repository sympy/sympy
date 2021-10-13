import pytest

from sympy import symbols
from sympy.tensor.array.expressions.array_expressions import ArrayElement, ArraySymbol


class TestArraySymbol:
    @pytest.mark.parametrize(
        "shape",
        [
            (),
            (3, 2, 4),
            symbols("k m n"),
        ],
    )
    def test_constructor(self, shape: tuple):
        A = ArraySymbol("A", *shape)
        assert A.name == "A"
        assert A.shape == shape

    def test_equality(self):
        m, n = symbols("m n")
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
        with pytest.raises(ValueError, match="shape is out of bounds"):
            A[4, 1, 3]
        with pytest.raises(ValueError, match="shape contains negative values"):
            A[0, -1]
        with pytest.raises(
            IndexError,
            match=(
                f"Too many indices for {ArrayElement.__name__}: parent "
                f"{ArraySymbol.__name__} is 3-dimensional, but 4 indices were given"
            )
        ):
            A[0, 1, 2, 5]
