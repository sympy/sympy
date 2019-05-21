from sympy.tensor.array.array_comprehension import ArrayComprehension
from sympy.abc import i, j, k


def test_array_comprehension():
    a = ArrayComprehension(i, (i, 1, 5))
    b = ArrayComprehension(i, (i, 1, j))
    assert a.doit() == [1, 2, 3, 4, 5]
    assert isinstance(b.doit(), ArrayComprehension)
    assert not isinstance(a.doit(), ArrayComprehension)