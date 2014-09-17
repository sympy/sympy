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

def test_reversed_rule_dispatch_sort_key():
    from sympy.core.compatibility import reversed_rule_dispatch_sort_key as key
    from sympy import symbols, Wild, S, sin, Function, cos

    x, y, z = symbols('x y z')
    u_, v_, w_ = map(Wild, list('uvw'))
    f = Function('f')

    pl = [x, u_, S.One, sin(x), sin(u_), sin(S.One)]
    res = [u_, x, S.One, sin(u_), sin(x), sin(S.One)]
    assert sorted(pl, key=key) == res

    pl = [f(x), f(u_), f(S.One)]
    res = [f(u_), f(x), f(S.One)]
    assert sorted(pl, key=key) == res

    pl = [f(x), f(u_), f(S.One), x, u_, S.One]
    res = [u_, x, S.One, f(u_), f(x), f(S.One)]
    assert sorted(pl, key=key) == res

    pl = [f(x), f(x, y), f(x, y, z)]
    res = [f(x), f(x, y), f(x, y, z)]
    assert sorted(pl, key=key) == res

    pl = [f(sin(x)),
         f(w_),
         f(sin(w_)),
         f(x, y),
         f(x, w_),
         f(v_, w_),
         f(x, y, z),
         f(x, y, v_),
         f(x, v_, y),
         f(x, sin(w_), cos(y))
        ]
    # Mathics output:
    #res = [
         #f(w_),
         #f(v_, w_),
         #f(sin(w_)),
         #f(sin(x)),
         #f(x, w_),
         #f(x, v_, y),
         #f(x, sin(w_), cos(y)),
         #f(x, y),
         #f(x, y, v_),
         #f(x, y, z),
    #]
    res = [
         f(w_),
         f(sin(w_)),
         f(sin(x)),
         f(v_, w_),
         f(x, w_),
         f(x, y),
         f(x, v_, y),
         f(x, y, v_),
         f(x, y, z),
         f(x, sin(w_), cos(y)),
    ]
    assert sorted(pl, key=key) == res

    pl = [f(x, y), f(u_, v_), f(sin(x), cos(y)), f(sin(u_), cos(v_))]
    res = [f(u_, v_), f(x, y), f(sin(u_), cos(v_)), f(sin(x), cos(y))]
    assert sorted(pl, key=key) == res
