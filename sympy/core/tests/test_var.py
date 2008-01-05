# Tests for var are in their own file, because var pollutes global namespace.

from sympy import Symbol, var
import py

# make z1 with call-depth = 1
def make_z1():
    var("z1")

# make z2 with call-depth = 2
def __make_z2():
    var("z2")

def make_z2():
    __make_z2()


def test_var():
    var("a")
    assert a  == Symbol("a")

    var("b bb cc zz _x")
    assert b  == Symbol("b")
    assert bb == Symbol("bb")
    assert cc == Symbol("cc")
    assert zz == Symbol("zz")
    assert _x == Symbol("_x")

    v = var(['d','e','fg'])
    assert d  == Symbol('d')
    assert e  == Symbol('e')
    assert fg == Symbol('fg')

    # check return value
    assert v  == (d, e, fg)

    # see if var() really injects into global namespace
    py.test.raises(NameError, "z1")
    make_z1()
    assert z1 == Symbol("z1")

    py.test.raises(NameError, "z2")
    make_z2()
    assert z2 == Symbol("z2")

def test_var_return():
    v1 = var('')
    v2 = var('q')
    v3 = var('q p')

    assert v1 == None
    assert v2 == Symbol('q')
    assert v3 == (Symbol('q'), Symbol('p'))
