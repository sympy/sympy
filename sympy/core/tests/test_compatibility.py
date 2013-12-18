from sympy.core.compatibility import default_sort_key, as_int, ordered
from sympy.core.singleton import S
from sympy.utilities.pytest import raises

from sympy.abc import x


def test_default_sort_key():
    func = lambda x: x
    assert sorted([func, x, func], key=default_sort_key) == [func, func, x]


def test_as_int():
    raises(ValueError, lambda : as_int(1.1))
    raises(ValueError, lambda : as_int([]))


def test_ordered():
    # Issue 4111 - this had been failing with python2/3 problems
    assert(list(ordered([{1:3, 2:4, 9:10}, {1:3}])) == \
               [{1: 3}, {1: 3, 2: 4, 9: 10}])
