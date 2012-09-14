import itertools

from sympy import (Add, Mul, Pow, Symbol, exp, sqrt, symbols, sympify, cse,
    Matrix, S, cos, sin, Eq)
from sympy.functions.special.hyper import meijerg
from sympy.simplify import cse_main, cse_opts
from sympy.utilities.pytest import XFAIL

w,x,y,z = symbols('w,x,y,z')
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
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]

def test_cse_single2():
    # Simple substitution, test for being able to pass the expression directly
    e = Add(Pow(x+y,2), sqrt(x+y))
    substs, reduced = cse(e, optimizations=[])
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]
    assert isinstance(cse(Matrix([[1]]))[1][0], Matrix)

def test_cse_not_possible():
    # No substitution possible.
    e = Add(x, y)
    substs, reduced = cse([e], optimizations=[])
    assert substs == []
    assert reduced == [x + y]
    # issue 3230
    eq = (meijerg((1, 2), (y, 4), (5,), [], x) + \
          meijerg((1, 3), (y, 4), (5,), [], x))
    assert cse(eq) == ([], [eq])

def test_nested_substitution():
    # Substitution within a substitution.
    e = Add(Pow(w*x+y,2), sqrt(w*x+y))
    substs, reduced = cse([e], optimizations=[])
    assert substs == [(x0, w*x+y)]
    assert reduced == [sqrt(x0) + x0**2]

def test_subtraction_opt():
    # Make sure subtraction is optimized.
    e = (x-y)*(z-y) + exp((x-y)*(z-y))
    substs, reduced = cse([e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, x - y), (x1, y - z), (x2, x0*x1)]
    assert reduced == [-x2 + exp(-x2)]
    assert cse(-(x - y)*(z - y) + exp(-(x - y)*(z - y))) == \
        ([(x0, (x - y)*(y - z))], [x0 + exp(x0)])
    # issue 978
    n = -1 + 1/x
    e = n/x/(-n)**2 - 1/n/x
    assert cse(e) == ([], [0])

def test_multiple_expressions():
    e1 = (x+y)*z
    e2 = (x+y)*w
    substs, reduced = cse([e1, e2], optimizations=[])
    assert substs == [(x0, x+y)]
    assert reduced == [x0*z, x0*w]
    l = [w*x*y + z, w*y]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == rsubsts
    assert reduced == [z + x*x0, x0]
    l = [w*x*y, w*x*y + z, w*y]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == rsubsts
    assert reduced == [x1, x1 + z, x0]
    l = [(x - z)*(y - z), x - z, y - z]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    substitutions = [
        [(x0, x - z), (x1, y - z)],
        [(x0, y - z), (x1, x - z)],
    ]
    assert substs in substitutions
    assert rsubsts in substitutions
    assert reduced == [x0*x1, x0, x1]
    l = [w*y + w + x + y + z, w*x*y]
    assert cse(l) == ([(x0, w*y)], [w + x + x0 + y + z, x*x0])
    assert cse([x + y, x + y + z]) == ([(x0, x + y)], [x0, z + x0])
    assert cse([x + y, x + z]) == ([], [x + y, x + z])
    assert cse([x*y, z + x*y , x*y*z + 3]) == \
           ([(x0, x*y)], [x0, z + x0, 3 + x0*z])
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*C]
    assert cse(l) == ([], l)
    l = [A*B*C, A*B]
    assert cse(l) == ([(x0, A*B)], [x0*C, x0])

@XFAIL
def test_powers():
    assert cse(x*y**2 + x*y) == ([(x0, x*y)], [x0*y + x0])

def test_issues_1399():
    assert cse(w/(x - y) + z/(y - x)) == ([], [(w - z)/(x - y)])

def test_issue_921():
    assert cse(x**5 + x**4 + x**3 + x**2) == ([(x0, x**2)], [x0*(x**3 + x + x0 + 1)])

def test_issue_1104():
    assert cse(sin(x**x)/x**x) == ([(x0, x**x)], [sin(x0)/x0])

def test_issue_3164():
    e = Eq(x*(-x + 1) + x*(x - 1), 0)
    assert cse(e) == ([], [True])

def test_dont_cse_tuples():
    from sympy import Subs, Function
    f = Function("f")
    g = Function("g")

    name_val, (expr,) = cse(
            Subs(f(x, y), (x, y), (0, 1))
            + Subs(g(x, y), (x, y), (0, 1)))

    assert name_val == []
    assert expr == (Subs(f(x, y), (x, y), (0, 1))
            + Subs(g(x, y), (x, y), (0, 1)))

    name_val, (expr,) = cse(
            Subs(f(x, y), (x, y), (0, x + y))
            + Subs(g(x, y), (x, y), (0, x + y)))

    assert name_val == [(x0, x + y)]
    assert expr == Subs(f(x, y), (x, y), (0, x0)) + \
          Subs(g(x, y), (x, y), (0, x0))

def test_pow_invpow():
    assert cse(1/x**2 + x**2) == \
        ([(x0, x**2)], [x0 + 1/x0])
    assert cse(x**2 + (1 + 1/x**2)/x**2) == \
        ([(x0, x**2)], [x0 + (1 + 1/x0)/x0])
    assert cse(1/x**2 + (1 + 1/x**2)*x**2) == \
        ([(x0, x**2)], [x0*(1 + 1/x0) + 1/x0])
    assert cse(cos(1/x**2) + sin(1/x**2)) == \
        ([(x0, x**2)], [sin(1/x0) + cos(1/x0)])
    assert cse(cos(x**2) + sin(x**2)) == \
        ([(x0, x**2)], [sin(x0) + cos(x0)])
    assert cse(y/(2 + x**2) + z/x**2/y) == \
        ([(x0, x**2)], [y/(x0 + 2) + z/(x0*y)])
    assert cse(exp(x**2) + x**2*cos(1/x**2)) == \
        ([(x0, x**2)], [x0*cos(1/x0) + exp(x0)])
    assert cse((1 + 1/x**2)/x**2) == \
        ([(x0, x**2)], [(1 + 1/x0)/x0])
    assert cse(x**(2*y) + x**(-2*y)) == \
        ([(x0, x**(2*y))], [x0 + 1/x0])

def test_postprocess():
    eq = (x + 1 + exp((x + 1)/(y + 1)) + cos(y + 1))
    assert cse([eq, Eq(x, z + 1), z - 2],
           postprocess=cse_main.cse_separate) == \
    [[(x0, y + 1), (x, z + 1), (x1, x + 1)], [x1 + exp(x1/x0) + cos(x0), z - 2]]
