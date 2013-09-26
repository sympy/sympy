from sympy.core.rules import Transform
from sympy import symbols, Wild, sin, cos
from sympy.core.rules import MapMatcher

from sympy.utilities.pytest import raises


def test_Transform():
    add1 = Transform(lambda x: x + 1, lambda x: x % 2 == 1)
    assert add1[1] == 2
    assert (1 in add1) is True
    assert add1.get(1) == 2

    raises(KeyError, lambda: add1[2])
    assert (2 in add1) is False
    assert add1.get(2) is None


def test_MapMatcher():
    x, y, z = symbols('x, y, z')
    a = Wild('a')
    b = Wild('b')
    c = Wild('c')

    mm = MapMatcher()
    mm[x+sin(a)] = lambda d: ("1", d)
    mm[x+y+a] = lambda d: ("2", d)
    mm[a/b] = lambda d: ("3", d)

    assert mm[x+sin(3)] == ("1", {a: 3})
    assert mm[x/y] == ("3", {a: x, b: y})
    assert mm[x+y+3] == ("2", {a: 3})
    assert mm[x+y] == ("2", {a: 0})
    assert mm[cos(x)+3] == ('3', {b: 1/(cos(x) + 3), a: 1})

    mm2 = MapMatcher()
    mm2[x+a] = lambda d : d
    mm2[2*x+2*a] = lambda d : d

    assert mm2[x+y] == {a: y}
    assert mm2[2*x+y] == {a: x+y}
    raises(KeyError, lambda: mm2[y])
