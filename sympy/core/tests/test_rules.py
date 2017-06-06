from sympy.core.rules import Transform
from sympy import symbols, Wild, sin, cos
from sympy.core.rules import MapMatcher
from sympy import Tuple
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
    mm[x+sin(a)] = Tuple(1, a, b)
    mm[x+y+a] = Tuple(2, a, b)
    mm[a/b] = Tuple(3, a, b)

    assert mm[x+sin(3)] == Tuple(1, 3, b)
    assert mm[x/y] == Tuple(3, x, y)
    assert mm[x+y+3] == Tuple(2, 3, b)
    assert mm[x+y] == Tuple(2, 0, b)
    assert mm[cos(x)+3] == Tuple(3, 1, 1/(cos(x) + 3))

    mm2 = MapMatcher()
    mm2[x+a] = Tuple(1, a)
    mm2[2*x+2*a] = Tuple(2, a)

    assert mm2[x+y] == Tuple(1, y)
    assert mm2[2*x+y] == Tuple(1, x+y)
    raises(KeyError, lambda: mm2[y])
