from sympy import Symbol, var, Function, FunctionClass
from sympy.utilities.pytest import raises

def test_var():
    ns = {"var": var, "raises": raises}
    eval("var('a')", ns)
    assert ns["a"] == Symbol("a")

    eval("var('b bb cc zz _x')", ns)
    assert ns["b"] == Symbol("b")
    assert ns["bb"] == Symbol("bb")
    assert ns["cc"] == Symbol("cc")
    assert ns["zz"] == Symbol("zz")
    assert ns["_x"] == Symbol("_x")

    v = eval("var(['d', 'e', 'fg'])", ns)
    assert ns['d'] == Symbol('d')
    assert ns['e'] == Symbol('e')
    assert ns['fg'] == Symbol('fg')

# check return value
    assert v == ['d', 'e', 'fg']

# make z1 with call-depth = 1

def _make_z1():
    eval("var('z1')", ns)

# make z2 with call-depth = 2

def __make_z2():
    eval("var('z2')", ns)

def _make_z2():
    __make_z2()

# see if var() really injects into global namespace
    "raises(NameError, lambda: z1)"
    _make_z1()
    assert ns["z1"] == Symbol("z1")

    "raises(NameError, lambda: z2)"
    _make_z2()
    assert ns["z2"] == Symbol("z2")



def test_var_return():
    ns = {"var": var, "raises": raises}
    "raises(ValueError, lambda: var(''))"
    v2 = eval("var('q')", ns)
    v3 = eval("var('q p')", ns)

    assert v2 == Symbol('q')
    assert v3 == (Symbol('q'), Symbol('p'))



def test_var_accepts_comma():
    ns = {"var": var}
    v1 = eval("var('x y z')", ns)
    v2 = eval("var('x,y,z')", ns)
    v3 = eval("var('x,y z')", ns)

    assert v1 == v2
    assert v1 == v3


def test_var_keywords():
    ns = {"var": var}
    eval("var('x y', real=True)", ns)
    assert ns['x'].is_real and ns['y'].is_real


def test_var_cls():
    ns = {"var": var, "isinstance": isinstance, "Function": Function, "FunctionClass": FunctionClass}
    f = eval("var('f', cls=Function)", ns)

    assert "isinstance(f, FunctionClass)"

    g, h = eval("var('g,h', cls=Function)", ns)

    assert "isinstance(g, FunctionClass)"
    assert "isinstance(h, FunctionClass)"
