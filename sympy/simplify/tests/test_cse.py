import itertools

from sympy import Add, Mul, Pow, Symbol, sin, sqrt, symbols, sympify, cse
from sympy.simplify import cse_main, cse_opts

w,x,y,z = symbols('wxyz')
x0,x1,x2 = list(itertools.islice(cse_main.numbered_symbols(), 0, 3))
negone = sympify(-1)


def test_numbered_symbols():
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(ns, 0, 10)) == [Symbol('y%s'%i) for i in range(0, 10)]
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(ns, 10, 20)) == [Symbol('y%s'%i) for i in range(10, 20)]
    ns = cse_main.numbered_symbols()
    assert list(itertools.islice(ns, 0, 10)) == [Symbol('x%s'%i) for i in range(0, 10)]

# Dummy "optimization" functions for testing.

def opt1(expr):
    return expr+y
def opt2(expr):
    return expr*z

def test_preprocess_for_cse():
    assert cse_main.preprocess_for_cse(x, [(opt1, None)]) == x+y
    assert cse_main.preprocess_for_cse(x, [(None, opt1)]) == x
    assert cse_main.preprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.preprocess_for_cse(x, [(opt1, opt2)]) == x+y
    assert cse_main.preprocess_for_cse(x, [(opt1, None), (opt2, None)]) == (x+y)*z

def test_postprocess_for_cse():
    assert cse_main.postprocess_for_cse(x, [(opt1, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(None, opt1)]) == x+y
    assert cse_main.postprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(opt1, opt2)]) == x*z
    # Note the reverse order of application.
    assert cse_main.postprocess_for_cse(x, [(None, opt1), (None, opt2)]) == x*z+y

def test_cse_single():
    # Simple substitution.
    e = Add(Pow(x+y,2), sqrt(x+y))
    substs, reduced = cse([e], optimizations=[])
    assert substs == [(x0, x+y)]
    assert reduced == [sqrt(x0) + x0**2]

def test_cse_single2():
    # Simple substitution, test for being able to pass the expression directly
    e = Add(Pow(x+y,2), sqrt(x+y))
    substs, reduced = cse(e, optimizations=[])
    assert substs == [(x0, x+y)]
    assert reduced == [sqrt(x0) + x0**2]

def test_cse_not_possible():
    # No substitution possible.
    e = Add(x,y)
    substs, reduced = cse([e], optimizations=[])
    assert substs == []
    assert reduced == [x+y]

def test_nested_substitution():
    # Substitution within a substitution.
    e = Add(Pow(w*x+y,2), sqrt(w*x+y))
    substs, reduced = cse([e], optimizations=[])
    assert substs == [(x0, w*x), (x1, x0+y)]
    assert reduced == [sqrt(x1) + x1**2]

def test_subtraction_opt():
    # Make sure subtraction is optimized.
    e = (x-y)*(z-y) + sin((x-y)*(z-y))
    substs, reduced = cse([e], optimizations=[(cse_opts.sub_pre,cse_opts.sub_post)])
    assert substs == [(x0, z-y), (x1, x-y), (x2, x0*x1)]
    assert reduced == [x2 + sin(x2)]

def test_multiple_expressions():
    e1 = (x+y)*z
    e2 = (x+y)*w
    substs, reduced = cse([e1, e2], optimizations=[])
    assert substs == [(x0, x+y)]
    assert reduced == [x0*z, x0*w]



