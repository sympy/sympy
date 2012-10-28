from sympy.core.compatibility import oset, default_sort_key
from sympy.core.singleton import S
from sympy.utilities.pytest import raises

from sympy.abc import x

def test_oset():
    a = oset([2, 1, 3, 1, x])
    b = oset([2, 1, 3, x])
    assert a == b
    assert len(b) == 4
    assert all(i in b for i in [1, 2, x])
    b.add(0)
    assert b == oset([2, 1, 3, x, 0])
    b.discard(1)
    assert b == oset([2, 3, x, 0])
    assert list(b) == [2, 3, x, 0]
    assert list(reversed(b)) == list(reversed([2, 3, x, 0]))
    z = b.pop()
    assert z == 0
    assert z is not S.Zero  # oset doesn't (and shouldn't) sympify
    assert b.pop(0) == 2
    assert repr(b) == 'oset([3, x])'
    assert repr(oset()) == 'oset()'
    assert b != [3, x]
    assert b == oset([3, x])
    assert a != b
    b.update([0, 1, 3])
    assert b == oset([3, x, 0, 1])

    a = oset([1, 2, 3])
    b = oset([3, 1])
    assert a.intersection(b) == oset([1, 3])
    assert b.intersection(a) == oset([3, 1])
    assert a.difference(b) == oset([2])
    assert a.difference([2]) == oset([1, 3])

    raises(KeyError, lambda: oset().pop(1))

def test_default_sort_key():
    func = lambda x: x
    assert sorted([func, x, func], key=default_sort_key) == [func, func, x]
