from sympy.tensor.array.array_comprehension import ArrayComprehension
from sympy.tensor.array import ImmutableDenseNDimArray
from sympy.abc import i, j, k, l
from sympy.utilities.pytest import raises


def test_array_comprehension():
    a = ArrayComprehension(i*j, (i, 1, 3), (j, 2, 4))
    b = ArrayComprehension(i, (i, 1, j))
    c = ArrayComprehension(i+j+k+l, (i, 1, 2), (j, 1, 3), (k, 1, 4), (l, 1, 5))
    assert a.doit().tolist() == [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
    assert isinstance(b.doit(), ArrayComprehension)
    assert isinstance(a.doit(), ImmutableDenseNDimArray)
    #assert b.subs(j, 3) == [3, 6, 9]
    assert c.doit().tolist() == [[[[4, 5, 6, 7, 8], [5, 6, 7, 8, 9], [6, 7, 8, 9, 10], [7, 8, 9, 10, 11]],
                                  [[5, 6, 7, 8, 9], [6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12]],
                                  [[6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12], [9, 10, 11, 12, 13]]],
                                 [[[5, 6, 7, 8, 9], [6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12]],
                                  [[6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12], [9, 10, 11, 12, 13]],
                                  [[7, 8, 9, 10, 11], [8, 9, 10, 11, 12], [9, 10, 11, 12, 13], [10, 11, 12, 13, 14]]]]
    raises(TypeError, lambda:ArrayComprehension(i*j, (i, 1, 3), (j, 2, [1, 3, 2])))
    raises(ValueError, lambda:ArrayComprehension(i*j, (i, 1, 3), (j, 2, 1)))
    raises(ValueError, lambda:ArrayComprehension(i*j, (i, 1, 3), (l, 1, 2)))
