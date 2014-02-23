import itertools

from sympy import (Add, Pow, Symbol, exp, sqrt, symbols, sympify, cse,
                   Matrix, S, cos, sin, Eq, Function, Tuple, RootOf, IndexedBase, Idx, MatrixSymbol)
from sympy.simplify.cse_opts import sub_pre, sub_post
from sympy.functions.special.hyper import meijerg
from sympy.simplify import cse_main, cse_opts
from sympy.utilities.pytest import XFAIL

w, x, y, z = symbols('w,x,y,z')
x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = symbols('x:13')


def test_numbered_symbols():
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(
        ns, 0, 10)) == [Symbol('y%s' % i) for i in range(0, 10)]
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(
        ns, 10, 20)) == [Symbol('y%s' % i) for i in range(10, 20)]
    ns = cse_main.numbered_symbols()
    assert list(itertools.islice(
        ns, 0, 10)) == [Symbol('x%s' % i) for i in range(0, 10)]

# Dummy "optimization" functions for testing.


def opt1(expr):
    return expr + y


def opt2(expr):
    return expr*z


def test_preprocess_for_cse():
    assert cse_main.preprocess_for_cse(x, [(opt1, None)]) == x + y
    assert cse_main.preprocess_for_cse(x, [(None, opt1)]) == x
    assert cse_main.preprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.preprocess_for_cse(x, [(opt1, opt2)]) == x + y
    assert cse_main.preprocess_for_cse(
        x, [(opt1, None), (opt2, None)]) == (x + y)*z


def test_postprocess_for_cse():
    assert cse_main.postprocess_for_cse(x, [(opt1, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(None, opt1)]) == x + y
    assert cse_main.postprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(opt1, opt2)]) == x*z
    # Note the reverse order of application.
    assert cse_main.postprocess_for_cse(
        x, [(None, opt1), (None, opt2)]) == x*z + y


def test_cse_single():
    # Simple substitution.
    e = Add(Pow(x + y, 2), sqrt(x + y))
    substs, reduced = cse([e])
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]


def test_cse_single2():
    # Simple substitution, test for being able to pass the expression directly
    e = Add(Pow(x + y, 2), sqrt(x + y))
    substs, reduced = cse(e)
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]
    assert isinstance(cse(Matrix([[1]]))[1][0], Matrix)


def test_cse_not_possible():
    # No substitution possible.
    e = Add(x, y)
    substs, reduced = cse([e])
    assert substs == []
    assert reduced == [x + y]
    # issue 3230
    eq = (meijerg((1, 2), (y, 4), (5,), [], x) +
          meijerg((1, 3), (y, 4), (5,), [], x))
    assert cse(eq) == ([], [eq])


def test_nested_substitution():
    # Substitution within a substitution.
    e = Add(Pow(w*x + y, 2), sqrt(w*x + y))
    substs, reduced = cse([e])
    assert substs == [(x0, w*x + y)]
    assert reduced == [sqrt(x0) + x0**2]


def test_subtraction_opt():
    # Make sure subtraction is optimized.
    e = (x - y)*(z - y) + exp((x - y)*(z - y))
    substs, reduced = cse(
        [e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, (x - y)*(y - z))]
    assert reduced == [-x0 + exp(-x0)]
    e = -(x - y)*(z - y) + exp(-(x - y)*(z - y))
    substs, reduced = cse(
        [e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, (x - y)*(y - z))]
    assert reduced == [x0 + exp(x0)]
    # issue 978
    n = -1 + 1/x
    e = n/x/(-n)**2 - 1/n/x
    assert cse(e, optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)]) == \
        ([], [0])


def test_multiple_expressions():
    e1 = (x + y)*z
    e2 = (x + y)*w
    substs, reduced = cse([e1, e2])
    assert substs == [(x0, x + y)]
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
    assert substs == [(x0, -z), (x1, x + x0), (x2, x0 + y)]
    assert rsubsts == [(x0, -z), (x1, x0 + y), (x2, x + x0)]
    assert reduced == [x1*x2, x1, x2]
    l = [w*y + w + x + y + z, w*x*y]
    assert cse(l) == ([(x0, w*y)], [w + x + x0 + y + z, x*x0])
    assert cse([x + y, x + y + z]) == ([(x0, x + y)], [x0, z + x0])
    assert cse([x + y, x + z]) == ([], [x + y, x + z])
    assert cse([x*y, z + x*y, x*y*z + 3]) == \
        ([(x0, x*y)], [x0, z + x0, 3 + x0*z])

@XFAIL # CSE of non-commutative Mul terms is disabled
def test_non_commutative_cse():
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*C]
    assert cse(l) == ([], l)
    l = [A*B*C, A*B]
    assert cse(l) == ([(x0, A*B)], [x0*C, x0])

# Test if CSE of non-commutative Mul terms is disabled
def test_bypass_non_commutatives():
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*C]
    assert cse(l) == ([], l)
    l = [A*B*C, A*B]
    assert cse(l) == ([], l)
    l = [B*C, A*B*C]
    assert cse(l) == ([], l)

@XFAIL # CSE fails when replacing non-commutative sub-expressions
def test_non_commutative_order():
    A, B, C = symbols('A B C', commutative=False)
    x0 = symbols('x0', commutative=False)
    l = [B+C, A*(B+C)]
    assert cse(l) == ([(x0, B+C)], [x0, A*x0])

@XFAIL
def test_powers():
    assert cse(x*y**2 + x*y) == ([(x0, x*y)], [x0*y + x0])


def test_issues_1399():
    assert cse(w/(x - y) + z/(y - x), optimizations='basic') == \
        ([], [(w - z)/(x - y)])


def test_issue_921():
    assert cse(x**5 + x**4 + x**3 + x**2, optimizations='basic') \
        == ([(x0, x**2)], [x0*(x**3 + x + x0 + 1)])


def test_issue_1104():
    assert cse(sin(x**x)/x**x) == ([(x0, x**x)], [sin(x0)/x0])


def test_issue_3164():
    e = Eq(x*(-x + 1) + x*(x - 1), 0)
    assert cse(e, optimizations='basic') == ([], [True])


def test_dont_cse_tuples():
    from sympy import Subs
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
        ([(x0, x**2), (x1, 1/x0)], [x0 + x1*(x1 + 1)])
    assert cse(1/x**2 + (1 + 1/x**2)*x**2) == \
        ([(x0, x**2), (x1, 1/x0)], [x0*(x1 + 1) + x1])
    assert cse(cos(1/x**2) + sin(1/x**2)) == \
        ([(x0, x**(-2))], [sin(x0) + cos(x0)])
    assert cse(cos(x**2) + sin(x**2)) == \
        ([(x0, x**2)], [sin(x0) + cos(x0)])
    assert cse(y/(2 + x**2) + z/x**2/y) == \
        ([(x0, x**2)], [y/(x0 + 2) + z/(x0*y)])
    assert cse(exp(x**2) + x**2*cos(1/x**2)) == \
        ([(x0, x**2)], [x0*cos(1/x0) + exp(x0)])
    assert cse((1 + 1/x**2)/x**2) == \
        ([(x0, x**(-2))], [x0*(x0 + 1)])
    assert cse(x**(2*y) + x**(-2*y)) == \
        ([(x0, x**(2*y))], [x0 + 1/x0])


def test_postprocess():
    eq = (x + 1 + exp((x + 1)/(y + 1)) + cos(y + 1))
    assert cse([eq, Eq(x, z + 1), z - 2, (z + 1)*(x + 1)],
        postprocess=cse_main.cse_separate) == \
        [[(x1, y + 1), (x2, z + 1), (x, x2), (x0, x + 1)],
        [x0 + exp(x0/x1) + cos(x1), z - 2, x0*x2]]


def test_issue1400():
    # previously, this gave 16 constants
    from sympy.abc import a, b
    B = Function('B')
    G = Function('G')
    t = Tuple(*
        (a, a + S(1)/2, 2*a, b, 2*a - b + 1, (sqrt(z)/2)**(-2*a + 1)*B(2*a -
        b, sqrt(z))*B(b - 1, sqrt(z))*G(b)*G(2*a - b + 1),
        sqrt(z)*(sqrt(z)/2)**(-2*a + 1)*B(b, sqrt(z))*B(2*a - b,
        sqrt(z))*G(b)*G(2*a - b + 1), sqrt(z)*(sqrt(z)/2)**(-2*a + 1)*B(b - 1,
        sqrt(z))*B(2*a - b + 1, sqrt(z))*G(b)*G(2*a - b + 1),
        (sqrt(z)/2)**(-2*a + 1)*B(b, sqrt(z))*B(2*a - b + 1,
        sqrt(z))*G(b)*G(2*a - b + 1), 1, 0, S(1)/2, z/2, -b + 1, -2*a + b,
        -2*a))
    c = cse(t)
    ans = (
        [(x0, 2*a), (x1, -b), (x2, x1 + 1), (x3, x0 + x2), (x4, sqrt(z)), (x5,
        B(x0 + x1, x4)), (x6, G(b)), (x7, G(x3)), (x8, -x0), (x9,
        (x4/2)**(x8 + 1)), (x10, x6*x7*x9*B(b - 1, x4)), (x11, x6*x7*x9*B(b,
        x4)), (x12, B(x3, x4))], [(a, a + S(1)/2, x0, b, x3, x10*x5,
        x11*x4*x5, x10*x12*x4, x11*x12, 1, 0, S(1)/2, z/2, x2, b + x8, x8)])
    assert ans == c


def test_issue_3070():
    r = RootOf(x**6 - 4*x**5 - 2, 1)
    assert cse(r) == ([], [r])
    # and a check that the right thing is done with the new
    # mechanism
    assert sub_post(sub_pre((-x - y)*z - x - y)) == -z*(x + y) - x - y

def test_cse_Indexed():
    len_y = 5
    y = IndexedBase('y', shape=(len_y,))
    x = IndexedBase('x', shape=(len_y,))
    Dy = IndexedBase('Dy', shape=(len_y-1,))
    i = Idx('i', len_y-1)

    expr1 = (y[i+1]-y[i])/(x[i+1]-x[i])
    expr2 = 1/(x[i+1]-x[i])
    replacements, reduced_exprs = cse([expr1, expr2])
    assert len(replacements) > 0

@XFAIL
def test_cse_MatrixSymbol():
    from sympy import MatrixSymbol
    A = MatrixSymbol('A', 3, 3)
    y = MatrixSymbol('y', 3, 1)

    expr1 = (A.T*A).I * A * y
    expr2 = (A.T*A) * A * y
    replacements, reduced_exprs = cse([expr1, expr2])
    assert len(replacements) > 0
