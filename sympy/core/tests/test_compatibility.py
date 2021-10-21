from sympy.core.compatibility import default_sort_key, as_int, ordered
from sympy.core.singleton import S
from sympy.testing.pytest import raises

from sympy.abc import x


def test_default_sort_key():
    func = lambda x: x
    assert sorted([func, x, func], key=default_sort_key) == [func, func, x]

    class C:
        def __repr__(self):
            return 'x.y'
    func = C()
    assert sorted([x, func], key=default_sort_key) == [func, x]


def test_as_int():
    raises(ValueError, lambda : as_int(True))
    raises(ValueError, lambda : as_int(1.1))
    raises(ValueError, lambda : as_int([]))
    raises(ValueError, lambda : as_int(S.NaN))
    raises(ValueError, lambda : as_int(S.Infinity))
    raises(ValueError, lambda : as_int(S.NegativeInfinity))
    raises(ValueError, lambda : as_int(S.ComplexInfinity))
    # for the following, limited precision makes int(arg) == arg
    # but the int value is not necessarily what a user might have
    # expected; Q.prime is more nuanced in its response for
    # expressions which might be complex representations of an
    # integer. This is not -- by design -- as_ints role.
    raises(ValueError, lambda : as_int(1e23))
    raises(ValueError, lambda : as_int(S('1.'+'0'*20+'1')))
    assert as_int(True, strict=False) == 1


def test_ordered():
    # Issue 7210 - this had been failing with python2/3 problems
    assert (list(ordered([{1:3, 2:4, 9:10}, {1:3}])) == \
               [{1: 3}, {1: 3, 2: 4, 9: 10}])
    # warnings should not be raised for identical items
    l = [1, 1]
    assert list(ordered(l, warn=True)) == l
    l = [[1], [2], [1]]
    assert list(ordered(l, warn=True)) == [[1], [1], [2]]
    raises(ValueError, lambda: list(ordered(['a', 'ab'], keys=[lambda x: x[0]],
        default=False, warn=True)))
