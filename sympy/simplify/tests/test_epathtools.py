"""Tests for tools for manipulation of expressions using paths. """

from sympy.simplify.epathtools import epath, EPath, eselect
from sympy.utilities.pytest import raises

from sympy import sin, cos, E
from sympy.abc import x, y, z, t

def test_eselect():
    expr = [((x, 1, t), 2), ((3, y, 4), z)]

    assert eselect(expr, "/*") == [((x, 1, t), 2), ((3, y, 4), z)]
    assert eselect(expr, "/*/*") == [(x, 1, t), 2, (3, y, 4), z]
    assert eselect(expr, "/*/*/*") == [x, 1, t, 3, y, 4]
    assert eselect(expr, "/*/*/*/*") == []

    assert eselect(expr, "/[:]") == [((x, 1, t), 2), ((3, y, 4), z)]
    assert eselect(expr, "/[:]/[:]") == [(x, 1, t), 2, (3, y, 4), z]
    assert eselect(expr, "/[:]/[:]/[:]") == [x, 1, t, 3, y, 4]
    assert eselect(expr, "/[:]/[:]/[:]/[:]") == []

    assert eselect(expr, "/*/[:]") == [(x, 1, t), 2, (3, y, 4), z]

    assert eselect(expr, "/*/[0]") == [(x, 1, t), (3, y, 4)]
    assert eselect(expr, "/*/[1]") == [2, z]
    assert eselect(expr, "/*/[2]") == []

    assert eselect(expr, "/*/int") == [2]
    assert eselect(expr, "/*/Symbol") == [z]
    assert eselect(expr, "/*/tuple") == [(x, 1, t), (3, y, 4)]
    assert eselect(expr, "/*/__iter__?") == [(x, 1, t), (3, y, 4)]

    assert eselect(expr, "/*/int|tuple") == [(x, 1, t), 2, (3, y, 4)]
    assert eselect(expr, "/*/Symbol|tuple") == [(x, 1, t), (3, y, 4), z]
    assert eselect(expr, "/*/int|Symbol|tuple") == [(x, 1, t), 2, (3, y, 4), z]

    assert eselect(expr, "/*/int|__iter__?") == [(x, 1, t), 2, (3, y, 4)]
    assert eselect(expr, "/*/Symbol|__iter__?") == [(x, 1, t), (3, y, 4), z]
    assert eselect(expr, "/*/int|Symbol|__iter__?") == [(x, 1, t), 2, (3, y, 4), z]

    assert eselect(expr, "/*/[0]/int") == [1, 3, 4]
    assert eselect(expr, "/*/[0]/Symbol") == [x, t, y]

    assert eselect(expr, "/*/[0]/int[1:]") == [1, 4]
    assert eselect(expr, "/*/[0]/Symbol[1:]") == [t, y]

    assert eselect(x + y + z + 1, "/Symbol") == [x, y, z]
    assert eselect(t + sin(x + 1) + cos(x + y + E), "/*/*/Symbol")   == [x, x, y]

def test_epath_apply():
    expr = [((x, 1, t), 2), ((3, y, 4), z)]
    func = lambda expr: expr**2

    assert epath("/*", expr, list) == [[(x, 1, t), 2], [(3, y, 4), z]]

    assert epath("/*/[0]", expr, list) == [([x, 1, t], 2), ([3, y, 4], z)]
    assert epath("/*/[1]", expr, func) == [((x, 1, t), 4), ((3, y, 4), z**2)]
    assert epath("/*/[2]", expr, list) == expr

    assert epath("/*/[0]/int", expr, func) == [((x, 1, t), 2), ((9, y, 16), z)]
    assert epath("/*/[0]/Symbol", expr, func) == [((x**2, 1, t**2), 2), ((3, y**2, 4), z)]
    assert epath("/*/[0]/int[1:]", expr, func) == [((x, 1, t), 2), ((3, y, 16), z)]
    assert epath("/*/[0]/Symbol[1:]", expr, func) == [((x, 1, t**2), 2), ((3, y**2, 4), z)]

    assert epath("/Symbol", x + y + z + 1, func) == x**2 + y**2 + z**2 + 1
    assert epath("/*/*/Symbol", t + sin(x + 1) + cos(x + y + E), func) == \
        t + sin(x**2 + 1) + cos(x**2 + y**2 + E)

def test_EPath():
    assert EPath("/*/[0]")._path == "/*/[0]"
    assert EPath(EPath("/*/[0]"))._path == "/*/[0]"
    assert isinstance(epath("/*/[0]"), EPath) is True

    assert repr(EPath("/*/[0]")) == "EPath('/*/[0]')"

    raises(ValueError, 'EPath("")')
    raises(ValueError, 'EPath("/")')
    raises(ValueError, 'EPath("/|x")')
    raises(ValueError, 'EPath("/[")')
    raises(ValueError, 'EPath("/[0]%")')

    raises(NotImplementedError, 'EPath("Symbol")')
