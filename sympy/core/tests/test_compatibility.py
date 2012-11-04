from sympy.core.compatibility import default_sort_key
from sympy.core.singleton import S
from sympy.utilities.pytest import raises

from sympy.abc import x

def test_default_sort_key():
    func = lambda x: x
    assert sorted([func, x, func], key=default_sort_key) == [func, func, x]
