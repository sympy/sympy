from sympy.core.compatibility import default_sort_key, as_int, ordered, iterable
from sympy.core.singleton import S
from sympy.utilities.pytest import raises

from sympy.abc import x


def test_default_sort_key():
    func = lambda x: x
    assert sorted([func, x, func], key=default_sort_key) == [func, func, x]


def test_as_int():
    raises(ValueError, lambda : as_int(1.1))
    raises(ValueError, lambda : as_int([]))


def test_iterable():
    assert iterable(0) is False
    assert iterable(1) is False
    assert iterable(None) is False


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
