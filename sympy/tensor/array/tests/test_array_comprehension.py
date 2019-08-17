from sympy.tensor.array.array_comprehension import ArrayComprehension
from sympy.tensor.array import ImmutableDenseNDimArray
from sympy.abc import i, j, k, l
from sympy.utilities.pytest import raises
from sympy.matrices import Matrix


def test_array_comprehension():
    a = ArrayComprehension(i*j, (i, 1, 3), (j, 2, 4))
    b = ArrayComprehension(i, (i, 1, j+1))
    c = ArrayComprehension(i+j+k+l, (i, 1, 2), (j, 1, 3), (k, 1, 4), (l, 1, 5))
    d = ArrayComprehension(k, (i, 1, 5))
    e = ArrayComprehension(i, (j, k+1, k+5))
    assert a.doit().tolist() == [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
    assert a.shape == (3, 3)
    assert a.is_shape_numeric == True
    assert a.tolist() == [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
    assert a.tomatrix() == Matrix([
                           [2, 3, 4],
                           [4, 6, 8],
                           [6, 9, 12]])
    assert len(a) == 9
    assert isinstance(b.doit(), ArrayComprehension)
    assert isinstance(a.doit(), ImmutableDenseNDimArray)
    assert b.subs(j, 3) == ArrayComprehension(i, (i, 1, 4))
    assert b.free_symbols == {j}
    assert b.shape == (j + 1,)
    assert b.rank() == 1
    assert b.is_shape_numeric == False
    assert c.free_symbols == set()
    assert c.function == i + j + k + l
    assert c.limits == ((i, 1, 2), (j, 1, 3), (k, 1, 4), (l, 1, 5))
    assert c.doit().tolist() == [[[[4, 5, 6, 7, 8], [5, 6, 7, 8, 9], [6, 7, 8, 9, 10], [7, 8, 9, 10, 11]],
                                  [[5, 6, 7, 8, 9], [6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12]],
                                  [[6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12], [9, 10, 11, 12, 13]]],
                                 [[[5, 6, 7, 8, 9], [6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12]],
                                  [[6, 7, 8, 9, 10], [7, 8, 9, 10, 11], [8, 9, 10, 11, 12], [9, 10, 11, 12, 13]],
                                  [[7, 8, 9, 10, 11], [8, 9, 10, 11, 12], [9, 10, 11, 12, 13], [10, 11, 12, 13, 14]]]]
    assert c.free_symbols == set()
    assert c.variables == [i, j, k, l]
    assert c.bound_symbols == [i, j, k, l]
    assert d.doit().tolist() == [k, k, k, k, k]
    assert len(e) == 5
    raises(TypeError, lambda: ArrayComprehension(i*j, (i, 1, 3), (j, 2, [1, 3, 2])))
    raises(ValueError, lambda: ArrayComprehension(i*j, (i, 1, 3), (j, 2, 1)))
    raises(ValueError, lambda: ArrayComprehension(i*j, (i, 1, 3), (j, 2, j+1)))
    raises(ValueError, lambda: len(ArrayComprehension(i*j, (i, 1, 3), (j, 2, j+4))))
    raises(TypeError, lambda: ArrayComprehension(i*j, (i, 0, i + 1.5), (j, 0, 2)))
    raises(ValueError, lambda: b.tolist())
    raises(ValueError, lambda: b.tomatrix())
    raises(ValueError, lambda: c.tomatrix())
