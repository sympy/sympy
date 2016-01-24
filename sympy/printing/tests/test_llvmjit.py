
from sympy.external import import_module


if import_module('llvmlite'):
    import sympy.printing.llvmjitcode as g
else:
    disabled = True

import sympy
from sympy.abc import a, b


# copied from numpy.isclose documentation
def isclose(a, b):
    rtol = 1e-5
    atol = 1e-8
    return abs(a-b) <= atol + rtol*abs(b)


def test_simple_expr():
    e = a + 1.0
    f = g.get_jit_callable(e, [a])
    res = float(e.subs({a: 4.0}).evalf())
    jit_res = f(4.0)

    assert isclose(jit_res, res)


def test_two_arg():
    e = 4.0*a + b + 3.0
    f = g.get_jit_callable(e, [a, b])
    res = float(e.subs({a: 4.0, b: 3.0}).evalf())
    jit_res = f(4.0, 3.0)

    assert isclose(jit_res, res)


def test_func():
    e = 4.0*sympy.exp(-a)
    f = g.get_jit_callable(e, [a])
    res = float(e.subs({a: 1.5}).evalf())
    jit_res = f(1.5)

    assert isclose(jit_res, res)
