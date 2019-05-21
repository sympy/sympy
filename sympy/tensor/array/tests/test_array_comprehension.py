from sympy.tensor.array.array_comprehension import ArrayComprehension
from sympy.abc import i, j


def test_array_comprehension():
    a = ArrayComprehension(i*j, (i, 1, 3), (j, 2, 4))
    b = ArrayComprehension(i, (i, 1, j))
    c = ArrayComprehension(i*j, (i, 1, 3))
    assert a.doit() == [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
    assert isinstance(b.doit(), ArrayComprehension)
    assert b.subs(j, 3) == [1, 2, 3]
    assert c.doit == [j, 2*j, 3*j]
    c.add_bound((j, 2, 4))
    assert c.doit() == [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
