# Tests for var are in their own file, because var pollutes global namespace.
import textwrap
from sympy import Symbol, var, Function, FunctionClass
from sympy.utilities.pytest import raises

# make z1 with call-depth = 1



# make z2 with call-depth = 2


def test_var():
    exec(textwrap.dedent("""
var("a")
assert a == Symbol("a")

var("b bb cc zz _x")
assert b == Symbol("b")
assert bb == Symbol("bb")
assert cc == Symbol("cc")
assert zz == Symbol("zz")
assert _x == Symbol("_x")

v = var(['d', 'e', 'fg'])
assert d == Symbol('d')
assert e == Symbol('e')
assert fg == Symbol('fg')

# check return value
assert v == [d, e, fg]

def _make_z1():
    z1 = var("z1")

def __make_z2():
    z2 = var("z2")

def _make_z2():
    __make_z2()

# see if var() really injects into global namespace
raises(NameError, lambda: z1)
_make_z1()
assert z1 == Symbol("z1")

raises(NameError, lambda: z2)
_make_z2()
assert z2 == Symbol("z2")
    """),
    {"var": var, "raises": raises, "Symbol": Symbol})


def test_var_return():
    exec("""
raises(ValueError, lambda: var(''))
v2 = var('q')
v3 = var('q p')

assert v2 == Symbol('q')
assert v3 == (Symbol('q'), Symbol('p'))
    """, {"var": var, "raises": raises, "Symbol": Symbol})


def test_var_accepts_comma():
    exec("""
v1 = var('x y z')
v2 = var('x,y,z')
v3 = var('x,y z')

assert v1 == v2
assert v1 == v3
    """, {"var": var})

def test_var_keywords():
    exec("""
var('x y', real=True)
assert x.is_real and y.is_real
    """, {"var": var})

def test_var_cls():
    exec("""
f = var('f', cls=Function)

assert isinstance(f, FunctionClass)

g, h = var('g,h', cls=Function)

assert isinstance(g, FunctionClass)
assert isinstance(h, FunctionClass)
    """, {"var": var, "isinstance": isinstance, "Function": Function, "FunctionClass": FunctionClass})
