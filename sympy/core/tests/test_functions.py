from sympy import Lambda, Symbol, Function, Derivative, sqrt, \
        log, exp, Rational, Float, sin, cos, acos, diff, I, re, im, \
        oo, zoo, nan, E, expand, pi, O, Sum, S, polygamma, loggamma
from sympy.utilities.pytest import XFAIL, raises
from sympy.abc import x, y, n
from sympy.core.function import PoleError
from sympy.utilities.iterables import subsets, variations
from sympy.core.compatibility import all

def test_f_expand_complex():
    f = Function('f')
    x = Symbol('x', real=True)
    z = Symbol('z')

    assert f(x).expand(complex=True)        == I*im(f(x)) + re(f(x))
    assert exp(x).expand(complex=True)      == exp(x)
    assert exp(I*x).expand(complex=True)    == cos(x) + I*sin(x)
    assert exp(z).expand(complex=True)      == cos(im(z))*exp(re(z)) + \
                                             I*sin(im(z))*exp(re(z))

def test_bug1():
    x = Symbol("x")
    w = Symbol("w")

    e = sqrt(-log(w))
    assert e.subs(log(w),-x) == sqrt(x)

    e = sqrt(-5*log(w))
    assert e.subs(log(w),-x) == sqrt(5*x)

def test_general_function():
    nu = Function('nu')
    x = Symbol("x")
    y = Symbol("y")

    e = nu(x)
    edx = e.diff(x)
    edy = e.diff(y)
    edxdx = e.diff(x).diff(x)
    edxdy = e.diff(x).diff(y)
    assert e == nu(x)
    assert edx != nu(x)
    assert edx == diff(nu(x), x)
    assert edy == 0
    assert edxdx == diff(diff(nu(x), x), x)
    assert edxdy == 0

def test_derivative_subs_bug():
    x = Symbol("x y")
    l = Function('l')
    n = Function('n')

    e = diff(n(x), x)
    assert e.subs(n(x), l(x)) != e
    assert e.subs(n(x), l(x)) == Derivative(l(x), x)
    assert e.subs(n(x), -l(x)) == Derivative(-l(x), x)

    assert e.subs(x, y) == Derivative(n(y), y)

def test_derivative_subs_self_bug():
    f = Function('f')
    d = diff(f(x), x)

    assert d.subs(d, y) == y


def test_derivative_linearity():
    x = Symbol("x")
    y = Symbol("y")
    n = Function('n')

    assert diff(-n(x), x) == -diff(n(x), x)
    assert diff(8*n(x), x) == 8*diff(n(x), x)
    assert diff(8*n(x), x) != 7*diff(n(x), x)
    assert diff(8*n(x)*x, x) == 8*n(x) + 8*x*diff(n(x), x)
    assert diff(8*n(x)*y*x, x) == 8*y*n(x) + 8*y*x*diff(n(x), x)

def test_derivative_evaluate():
    x = Symbol('x')

    assert Derivative(sin(x), x) != diff(sin(x), x)
    assert Derivative(sin(x), x).doit() == diff(sin(x), x)

    f = Function('f')
    assert Derivative(Derivative(f(x), x), x) == diff(f(x), x, x)
    assert Derivative(sin(x), x, 0) == sin(x)

def test_diff_symbols():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    f = Function('f')

    assert diff(f(x, y, z), x, y, z) == Derivative(f(x, y, z), x, y, z)
    assert diff(f(x, y, z), x, x, x) == Derivative(f(x, y, z), x, x, x)
    assert diff(f(x, y, z), x, 3) == Derivative(f(x, y, z), x, 3)

    # issue 1929
    assert [diff(-z + x/y, sym) for sym in (z, x, y)] == [-1, 1/y, -x/y**2]
    assert diff(f(x, y, z), x, y, z, 2) == Derivative(f(x, y, z), x, y, z, z)
    assert diff(f(x, y, z), x, y, z, 2, evaluate=False) == \
        Derivative(f(x, y, z), x, y, z, z)
    assert Derivative(f(x, y, z), x, y, z)._eval_derivative(z) == \
        Derivative(f(x, y, z), x, y, z, z)
    assert Derivative(Derivative(f(x, y, z), x), y)._eval_derivative(z) == \
        Derivative(f(x, y, z), x, y, z)

def test_Lambda():
    e = Lambda(x, x**2)
    f = Function('f')
    assert e(4) == 16
    assert e(x) == x**2
    assert e(y) == y**2

    assert Lambda(x, x**2) == Lambda(x, x**2)
    assert Lambda(x, x**2) == Lambda(y, y**2)
    assert Lambda(x, x**2) != Lambda(y, y**2+1)
    assert Lambda((x, y), x**y) == Lambda((y, x), y**x)
    assert Lambda((x, y), x**y) != Lambda((x, y), y**x)

    assert Lambda((x, y), x**y)(x, y) == x**y
    assert Lambda((x, y), x**y)(3, 3) == 3**3
    assert Lambda((x, y), x**y)(x, 3) == x**3
    assert Lambda((x, y), x**y)(3, y) == 3**y
    assert Lambda(x, f(x))(x) == f(x)
    assert Lambda(x, x**2)(e(x)) == x**4
    assert e(e(x)) == x**4

    assert Lambda((x, y), x+y).nargs == 2

    z = Symbol('z')
    t = Symbol('t')
    p = x, y, z, t
    assert Lambda(p, t*(x+y+z))(*p) == t * (x + y + z)

    assert Lambda(x, 2*x) + Lambda(y, 2*y) == 2*Lambda(x, 2*x)
    assert Lambda(x, 2*x) not in [ Lambda(x, x) ]

def test_IdentityFunction():
    assert Lambda(x, x) is Lambda(y, y) is S.IdentityFunction
    assert Lambda(x, 2*x) is not S.IdentityFunction
    assert Lambda((x, y), x) is not S.IdentityFunction

def test_Lambda_symbols():
    assert Lambda(x, 2*x).free_symbols == set()
    assert Lambda(x, x*y).free_symbols == set([y])

def test_Lambda_arguments():
    raises(TypeError, 'Lambda(x, 2*x)(x, y)')
    raises(TypeError, 'Lambda((x, y), x+y)(x)')

def test_Lambda_equality():
    assert Lambda(x, 2*x) != Lambda((x,y), 2*x)
    assert (Lambda(x, 2*x) == Lambda((x,y), 2*x)) is False
    assert Lambda((x, y), 2*x) == Lambda((x, y), 2*x)
    assert (Lambda((x, y), 2*x) != Lambda((x, y), 2*x)) is False
    assert Lambda(x, 2*x) != 2*x
    assert (Lambda(x, 2*x) == 2*x) is False


def test_expand_function():
    assert expand(x+y) == x + y
    assert expand(x+y, complex=True) == I*im(x) + I*im(y) + re(x) + re(y)
    assert expand((x + y)**11, modulus=11) == x**11 + y**11

def test_function_comparable():
    x = Symbol('x')

    assert sin(x).is_comparable == False
    assert cos(x).is_comparable == False

    assert sin(Float('0.1')).is_comparable   == True
    assert cos(Float('0.1')).is_comparable   == True

    assert sin(E).is_comparable     == True
    assert cos(E).is_comparable     == True

    assert sin(Rational(1,3)).is_comparable == True
    assert cos(Rational(1,3)).is_comparable == True

@XFAIL
def test_function_comparable():
    assert sin(oo).is_comparable    == False
    assert sin(-oo).is_comparable   == False
    assert sin(zoo).is_comparable   == False
    assert sin(nan).is_comparable   == False

@XFAIL
def test_deriv1():
    # This requires derivatives evaluated at a point.  This is written in such
    # a way that it will XPASS when that is implemented.  When it is, fix this
    # and test_deriv2() below.
    f=Function('f')
    g=Function('g')
    x = Symbol('x')
    f(g(x)).diff(x)

@XFAIL
def test_deriv2():
    # These all requre derivatives evaluated at a point (issue 1620) to work.
    # See issue 1525
    f = Function('f')
    x = Symbol('x')

    assert f(2*x).diff(x) == 2*Derivative(f(2*x), 2*x)
    assert (f(x)**3).diff(x) == 3*f(x)**2*f(x).diff(x)
    assert (f(2*x)**3).diff(x) == 6*f(2*x)**2*Derivative(f(2*x), 2*x)

    assert f(2+x).diff(x) == Derivative(f(2+x), 2+x)
    assert f(2+3*x).diff(x) == 3*Derivative(f(2+3*x), 2+3*x)
    assert f(sin(x)).diff(x) == Derivative(f(sin(x)), sin(x)) * cos(x)
    assert f(3*sin(x)).diff(x) == 3*Derivative(f(3*sin(x)), 3*sin(x)) * cos(x)

def test_deriv3():
    x = Symbol('x')

    assert (x**3).diff(x) == 3*x**2
    assert (x**3).diff(x, evaluate=False) != 3*x**2
    assert (x**3).diff(x, evaluate=False) == Derivative(x**3, x)

    assert diff(x**3, x) == 3*x**2
    assert diff(x**3, x, evaluate=False) != 3*x**2
    assert diff(x**3, x, evaluate=False) == Derivative(x**3, x)

def test_func_deriv():
    f = Function('f')
    x = Symbol('x')

    assert f(x).diff(x) == Derivative(f(x), x)

def test_suppressed_evaluation():
    a = sin(0, evaluate=False)
    assert a != 0
    assert a.func is sin
    assert a.args == (0,)

def test_function_evalf():
    def eq(a,b,eps):
        return abs(a-b) < eps
    assert eq(sin(1).evalf(15), Float("0.841470984807897"), 1e-13)
    assert eq(sin(2).evalf(25), Float("0.9092974268256816953960199",25), 1e-23)
    assert eq(sin(1+I).evalf(15), Float("1.29845758141598") + Float("0.634963914784736")*I, 1e-13)
    assert eq(exp(1+I).evalf(15), Float("1.46869393991588") + Float("2.28735528717884239")*I, 1e-13)
    assert eq(exp(-0.5+1.5*I).evalf(15), Float("0.0429042815937374") + Float("0.605011292285002")*I, 1e-13)
    assert eq(log(pi+sqrt(2)*I).evalf(15), Float("1.23699044022052") + Float("0.422985442737893")*I, 1e-13)
    assert eq(cos(100).evalf(15), Float("0.86231887228768"), 1e-13)

def test_extensibility_eval():
    class MyFunc(Function):
        @classmethod
        def eval(cls, *args):
            return (0,0,0)
    assert MyFunc(0) == (0,0,0)

def test_function_non_commutative():
    x = Symbol('x', commutative=False)
    f = Function('f')
    assert f(x).is_commutative == False
    assert sin(x).is_commutative == False
    assert exp(x).is_commutative == False
    assert log(x).is_commutative == False

def test_function__eval_nseries():
    x = Symbol('x')
    assert sin(x)._eval_nseries(x,2,None) == x + O(x**2)
    assert sin(x+1)._eval_nseries(x,2,None) == x*cos(1) + sin(1) + O(x**2)
    assert sin(pi*(1-x))._eval_nseries(x,2,None) == pi*x + O(x**2)
    assert acos(1-x**2)._eval_nseries(x,2,None) == sqrt(2)*x + O(x**2)
    assert polygamma(n,x+1)._eval_nseries(x,2,None) == \
                   polygamma(n,1) + polygamma(n+1,1)*x + O(x**2)
    raises(PoleError, 'sin(1/x)._eval_nseries(x,2,None)')
    raises(PoleError, 'acos(1-x)._eval_nseries(x,2,None)')
    raises(PoleError, 'acos(1+x)._eval_nseries(x,2,None)')
    assert loggamma(1/x)._eval_nseries(x,0,None) \
           == log(x)/2 - log(x)/x - 1/x + O(1, x)
    l = Symbol('l')
    assert loggamma(log(1/x)).nseries(x,n=1,logx=y) == loggamma(-y)

def test_doit():
    n = Symbol('n', integer = True)
    f = Sum(2 * n * x, (n, 1, 3))
    d = Derivative(f, x)
    assert d.doit() == 12
    assert d.doit(deep = False) == Sum(2*n, (n, 1, 3))

def test_evalf_default():
    from sympy.functions.special.gamma_functions import polygamma
    assert type(sin(4.0)) == Float
    assert type(re(sin(I + 1.0))) == Float
    assert type(im(sin(I + 1.0))) == Float
    assert type(sin(4)) == sin
    assert type(polygamma(2,4.0)) == Float
    assert type(sin(Rational(1,4))) == sin

def test_issue2300():
    args = [x, y, S(2), S.Half]
    def ok(a):
        """Return True if the input args for diff are ok"""
        if not a: return False
        if a[0].is_Symbol is False: return False
        s_at = [i for i in range(len(a)) if a[i].is_Symbol]
        n_at = [i for i in range(len(a)) if not a[i].is_Symbol]
        # every symbol is followed by symbol or int
        # every number is followed by a symbol
        return (all([a[i+1].is_Symbol or a[i+1].is_Integer
            for i in s_at if i+1<len(a)]) and
            all([a[i+1].is_Symbol
            for i in n_at if i+1<len(a)]))
    eq = x**10*y**8
    for a in subsets(args):
        for v in variations(a, len(a)):
            if ok(v):
                noraise = eq.diff(*v)
            else:
                raises(ValueError, 'eq.diff(*v)')
