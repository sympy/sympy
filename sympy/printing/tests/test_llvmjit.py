
from sympy.external import import_module
import ctypes


if import_module('llvmlite'):
    import sympy.printing.llvmjitcode as g
else:
    disabled = True

import sympy
from sympy.abc import a, b, n


# copied from numpy.isclose documentation
def isclose(a, b):
    rtol = 1e-5
    atol = 1e-8
    return abs(a-b) <= atol + rtol*abs(b)


def test_simple_expr():
    e = a + 1.0
    f = g.llvm_callable([a], e)
    res = float(e.subs({a: 4.0}).evalf())
    jit_res = f(4.0)

    assert isclose(jit_res, res)


def test_two_arg():
    e = 4.0*a + b + 3.0
    f = g.llvm_callable([a, b], e)
    res = float(e.subs({a: 4.0, b: 3.0}).evalf())
    jit_res = f(4.0, 3.0)

    assert isclose(jit_res, res)


def test_func():
    e = 4.0*sympy.exp(-a)
    f = g.llvm_callable([a], e)
    res = float(e.subs({a: 1.5}).evalf())
    jit_res = f(1.5)

    assert isclose(jit_res, res)


def test_two_func():
    e = 4.0*sympy.exp(-a) + sympy.exp(b)
    f = g.llvm_callable([a, b], e)
    res = float(e.subs({a: 1.5, b: 2.0}).evalf())
    jit_res = f(1.5, 2.0)

    assert isclose(jit_res, res)


def test_callback():
    e = a + 1.2
    f = g.llvm_callable([a], e, callback_type='scipy.integrate.test')
    m = ctypes.c_int(1)
    array_type = ctypes.c_double * 1
    inp = {a: 2.2}
    array = array_type(inp[a])
    jit_res = f(m, array)

    res = float(e.subs(inp).evalf())

    assert isclose(jit_res, res)


def test_callback_cubature():
    e = a + 1.2
    f = g.llvm_callable([a], e, callback_type='cubature')
    m = ctypes.c_int(1)
    array_type = ctypes.c_double * 1
    inp = {a: 2.2}
    array = array_type(inp[a])
    out_array = array_type(0.0)
    jit_ret = f(m, array, None, m, out_array)

    assert jit_ret == 0

    res = float(e.subs(inp).evalf())

    assert isclose(out_array[0], res)


def test_callback_two():
    e = 3*a*b
    f = g.llvm_callable([a, b], e, callback_type='scipy.integrate.test')
    m = ctypes.c_int(2)
    array_type = ctypes.c_double * 2
    inp = {a: 0.2, b: 1.7}
    array = array_type(inp[a], inp[b])
    jit_res = f(m, array)

    res = float(e.subs(inp).evalf())

    assert isclose(jit_res, res)


def test_callback_alt_two():
    d = sympy.IndexedBase('d')
    e = 3*d[0]*d[1]
    f = g.llvm_callable([n, d], e, callback_type='scipy.integrate.test')
    m = ctypes.c_int(2)
    array_type = ctypes.c_double * 2
    inp = {d[0]: 0.2, d[1]: 1.7}
    array = array_type(inp[d[0]], inp[d[1]])
    jit_res = f(m, array)

    res = float(e.subs(inp).evalf())

    assert isclose(jit_res, res)
