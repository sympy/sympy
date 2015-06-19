from sympy import (Lambda, Symbol, Function, Derivative, Subs, sqrt,
        log, exp, Rational, Float, sin, cos, acos, diff, I, re, im,
        E, expand, pi, O, Sum, S, polygamma, loggamma, expint,
        Tuple, Dummy, Eq, Expr, symbols, nfloat)
from sympy.utilities.pytest import XFAIL, raises
from sympy.abc import t, w, x, y, z
from sympy.core.function import PoleError, _mexpand
from sympy.sets.sets import FiniteSet
from sympy.solvers import solve
from sympy.utilities.iterables import subsets, variations
from sympy.core.cache import clear_cache
from sympy.core.compatibility import range

f, g, h = symbols('f g h', cls=Function)


def test_f_expand_complex():
    x = Symbol('x', real=True)

    assert f(x).expand(complex=True) == I*im(f(x)) + re(f(x))
    assert exp(x).expand(complex=True) == exp(x)
    assert exp(I*x).expand(complex=True) == cos(x) + I*sin(x)
    assert exp(z).expand(complex=True) == cos(im(z))*exp(re(z)) + \
        I*sin(im(z))*exp(re(z))


def test_bug1():
    e = sqrt(-log(w))
    assert e.subs(log(w), -x) == sqrt(x)

    e = sqrt(-5*log(w))
    assert e.subs(log(w), -x) == sqrt(5*x)


def test_general_function():
    nu = Function('nu')

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
    e = diff(g(x), x)
    assert e.subs(g(x), f(x)) != e
    assert e.subs(g(x), f(x)) == Derivative(f(x), x)
    assert e.subs(g(x), -f(x)) == Derivative(-f(x), x)

    assert e.subs(x, y) == Derivative(g(y), y)


def test_derivative_subs_self_bug():
    d = diff(f(x), x)

    assert d.subs(d, y) == y


def test_derivative_linearity():
    assert diff(-f(x), x) == -diff(f(x), x)
    assert diff(8*f(x), x) == 8*diff(f(x), x)
    assert diff(8*f(x), x) != 7*diff(f(x), x)
    assert diff(8*f(x)*x, x) == 8*f(x) + 8*x*diff(f(x), x)
    assert diff(8*f(x)*y*x, x) == 8*y*f(x) + 8*y*x*diff(f(x), x)


def test_derivative_evaluate():
    assert Derivative(sin(x), x) != diff(sin(x), x)
    assert Derivative(sin(x), x).doit() == diff(sin(x), x)

    assert Derivative(Derivative(f(x), x), x) == diff(f(x), x, x)
    assert Derivative(sin(x), x, 0) == sin(x)


def test_diff_symbols():
    assert diff(f(x, y, z), x, y, z) == Derivative(f(x, y, z), x, y, z)
    assert diff(f(x, y, z), x, x, x) == Derivative(f(x, y, z), x, x, x)
    assert diff(f(x, y, z), x, 3) == Derivative(f(x, y, z), x, 3)

    # issue 5028
    assert [diff(-z + x/y, sym) for sym in (z, x, y)] == [-1, 1/y, -x/y**2]
    assert diff(f(x, y, z), x, y, z, 2) == Derivative(f(x, y, z), x, y, z, z)
    assert diff(f(x, y, z), x, y, z, 2, evaluate=False) == \
        Derivative(f(x, y, z), x, y, z, z)
    assert Derivative(f(x, y, z), x, y, z)._eval_derivative(z) == \
        Derivative(f(x, y, z), x, y, z, z)
    assert Derivative(Derivative(f(x, y, z), x), y)._eval_derivative(z) == \
        Derivative(f(x, y, z), x, y, z)


def test_Function():
    class myfunc(Function):
        @classmethod
        def eval(cls, x):
            return

    assert myfunc.nargs == FiniteSet(1)
    assert myfunc(x).nargs == FiniteSet(1)
    raises(TypeError, lambda: myfunc(x, y).nargs)

    class myfunc(Function):
        @classmethod
        def eval(cls, *x):
            return

    assert myfunc.nargs == S.Naturals0
    assert myfunc(x).nargs == S.Naturals0

def test_nargs():
    f = Function('f')
    assert f.nargs == S.Naturals0
    assert f(1).nargs == S.Naturals0
    assert Function('f', nargs=2)(1, 2).nargs == FiniteSet(2)
    assert sin.nargs == FiniteSet(1)
    assert sin(2).nargs == FiniteSet(1)
    assert log.nargs == FiniteSet(1, 2)
    assert log(2).nargs == FiniteSet(1, 2)
    assert Function('f', nargs=2).nargs == FiniteSet(2)
    assert Function('f', nargs=0).nargs == FiniteSet(0)


def test_Lambda():
    e = Lambda(x, x**2)
    assert e(4) == 16
    assert e(x) == x**2
    assert e(y) == y**2

    assert Lambda(x, x**2) == Lambda(x, x**2)
    assert Lambda(x, x**2) == Lambda(y, y**2)
    assert Lambda(x, x**2) != Lambda(y, y**2 + 1)
    assert Lambda((x, y), x**y) == Lambda((y, x), y**x)
    assert Lambda((x, y), x**y) != Lambda((x, y), y**x)

    assert Lambda((x, y), x**y)(x, y) == x**y
    assert Lambda((x, y), x**y)(3, 3) == 3**3
    assert Lambda((x, y), x**y)(x, 3) == x**3
    assert Lambda((x, y), x**y)(3, y) == 3**y
    assert Lambda(x, f(x))(x) == f(x)
    assert Lambda(x, x**2)(e(x)) == x**4
    assert e(e(x)) == x**4

    assert Lambda((x, y), x + y).nargs == FiniteSet(2)

    p = x, y, z, t
    assert Lambda(p, t*(x + y + z))(*p) == t * (x + y + z)

    assert Lambda(x, 2*x) + Lambda(y, 2*y) == 2*Lambda(x, 2*x)
    assert Lambda(x, 2*x) not in [ Lambda(x, x) ]
    raises(ValueError, lambda: Lambda(1, x))
    assert Lambda(x, 1)(1) is S.One


def test_IdentityFunction():
    assert Lambda(x, x) is Lambda(y, y) is S.IdentityFunction
    assert Lambda(x, 2*x) is not S.IdentityFunction
    assert Lambda((x, y), x) is not S.IdentityFunction


def test_Lambda_symbols():
    assert Lambda(x, 2*x).free_symbols == set()
    assert Lambda(x, x*y).free_symbols == set([y])


def test_Lambda_arguments():
    raises(TypeError, lambda: Lambda(x, 2*x)(x, y))
    raises(TypeError, lambda: Lambda((x, y), x + y)(x))


def test_Lambda_equality():
    assert Lambda(x, 2*x) == Lambda(y, 2*y)
    # although variables are casts as Dummies, the expressions
    # should still compare equal
    assert Lambda((x, y), 2*x) == Lambda((x, y), 2*x)
    assert Lambda(x, 2*x) != Lambda((x, y), 2*x)
    assert Lambda(x, 2*x) != 2*x


def test_Subs():
    assert Subs(x, x, 0) == Subs(y, y, 0)
    assert Subs(x, x, 0).subs(x, 1) == Subs(x, x, 0)
    assert Subs(y, x, 0).subs(y, 1) == Subs(1, x, 0)
    assert Subs(f(x), x, 0).doit() == f(0)
    assert Subs(f(x**2), x**2, 0).doit() == f(0)
    assert Subs(f(x, y, z), (x, y, z), (0, 1, 1)) != \
        Subs(f(x, y, z), (x, y, z), (0, 0, 1))
    assert Subs(f(x, y), (x, y, z), (0, 1, 1)) == \
        Subs(f(x, y), (x, y, z), (0, 1, 2))
    assert Subs(f(x, y), (x, y, z), (0, 1, 1)) != \
        Subs(f(x, y) + z, (x, y, z), (0, 1, 0))
    assert Subs(f(x, y), (x, y), (0, 1)).doit() == f(0, 1)
    assert Subs(Subs(f(x, y), x, 0), y, 1).doit() == f(0, 1)
    raises(ValueError, lambda: Subs(f(x, y), (x, y), (0, 0, 1)))
    raises(ValueError, lambda: Subs(f(x, y), (x, x, y), (0, 0, 1)))

    assert len(Subs(f(x, y), (x, y), (0, 1)).variables) == 2
    assert Subs(f(x, y), (x, y), (0, 1)).point == Tuple(0, 1)

    assert Subs(f(x), x, 0) == Subs(f(y), y, 0)
    assert Subs(f(x, y), (x, y), (0, 1)) == Subs(f(x, y), (y, x), (1, 0))
    assert Subs(f(x)*y, (x, y), (0, 1)) == Subs(f(y)*x, (y, x), (0, 1))
    assert Subs(f(x)*y, (x, y), (1, 1)) == Subs(f(y)*x, (x, y), (1, 1))

    assert Subs(f(x), x, 0).subs(x, 1).doit() == f(0)
    assert Subs(f(x), x, y).subs(y, 0) == Subs(f(x), x, 0)
    assert Subs(y*f(x), x, y).subs(y, 2) == Subs(2*f(x), x, 2)
    assert (2 * Subs(f(x), x, 0)).subs(Subs(f(x), x, 0), y) == 2*y

    assert Subs(f(x), x, 0).free_symbols == set([])
    assert Subs(f(x, y), x, z).free_symbols == set([y, z])

    assert Subs(f(x).diff(x), x, 0).doit(), Subs(f(x).diff(x), x, 0)
    assert Subs(1 + f(x).diff(x), x, 0).doit(), 1 + Subs(f(x).diff(x), x, 0)
    assert Subs(y*f(x, y).diff(x), (x, y), (0, 2)).doit() == \
        2*Subs(Derivative(f(x, 2), x), x, 0)
    assert Subs(y**2*f(x), x, 0).diff(y) == 2*y*f(0)

    e = Subs(y**2*f(x), x, y)
    assert e.diff(y) == e.doit().diff(y) == y**2*Derivative(f(y), y) + 2*y*f(y)

    assert Subs(f(x), x, 0) + Subs(f(x), x, 0) == 2*Subs(f(x), x, 0)
    e1 = Subs(z*f(x), x, 1)
    e2 = Subs(z*f(y), y, 1)
    assert e1 + e2 == 2*e1
    assert e1.__hash__() == e2.__hash__()
    assert Subs(z*f(x + 1), x, 1) not in [ e1, e2 ]
    assert Derivative(
        f(x), x).subs(x, g(x)) == Subs(Derivative(f(x), x), (x,), (g(x),))
    assert Subs(f(x)*cos(y) + z, (x, y), (0, pi/3)).n(2) == \
        Subs(f(x)*cos(y) + z, (x, y), (0, pi/3)).evalf(2) == \
        z + Rational('1/2').n(2)*f(0)

    assert f(x).diff(x).subs(x, 0).subs(x, y) == f(x).diff(x).subs(x, 0)
    assert (x*f(x).diff(x).subs(x, 0)).subs(x, y) == y*f(x).diff(x).subs(x, 0)


@XFAIL
def test_Subs2():
    # this reflects a limitation of subs(), probably won't fix
    assert Subs(f(x), x**2, x).doit() == f(sqrt(x))


def test_expand_function():
    assert expand(x + y) == x + y
    assert expand(x + y, complex=True) == I*im(x) + I*im(y) + re(x) + re(y)
    assert expand((x + y)**11, modulus=11) == x**11 + y**11


def test_function_comparable():
    assert sin(x).is_comparable is False
    assert cos(x).is_comparable is False

    assert sin(Float('0.1')).is_comparable is True
    assert cos(Float('0.1')).is_comparable is True

    assert sin(E).is_comparable is True
    assert cos(E).is_comparable is True

    assert sin(Rational(1, 3)).is_comparable is True
    assert cos(Rational(1, 3)).is_comparable is True


@XFAIL
def test_function_comparable_infinities():
    assert sin(oo).is_comparable is False
    assert sin(-oo).is_comparable is False
    assert sin(zoo).is_comparable is False
    assert sin(nan).is_comparable is False


def test_deriv1():
    # These all requre derivatives evaluated at a point (issue 4719) to work.
    # See issue 4624
    assert f(2*x).diff(x) == 2*Subs(Derivative(f(x), x), Tuple(x), Tuple(2*x))
    assert (f(x)**3).diff(x) == 3*f(x)**2*f(x).diff(x)
    assert (
        f(2*x)**3).diff(x) == 6*f(2*x)**2*Subs(Derivative(f(x), x), Tuple(x),
            Tuple(2*x))

    assert f(2 + x).diff(x) == Subs(Derivative(f(x), x), Tuple(x), Tuple(x + 2))
    assert f(2 + 3*x).diff(x) == 3*Subs(Derivative(f(x), x), Tuple(x),
            Tuple(3*x + 2))
    assert f(3*sin(x)).diff(x) == 3*cos(x)*Subs(Derivative(f(x), x),
            Tuple(x), Tuple(3*sin(x)))


def test_deriv2():
    assert (x**3).diff(x) == 3*x**2
    assert (x**3).diff(x, evaluate=False) != 3*x**2
    assert (x**3).diff(x, evaluate=False) == Derivative(x**3, x)

    assert diff(x**3, x) == 3*x**2
    assert diff(x**3, x, evaluate=False) != 3*x**2
    assert diff(x**3, x, evaluate=False) == Derivative(x**3, x)


def test_func_deriv():
    assert f(x).diff(x) == Derivative(f(x), x)
    # issue 4534
    assert f(x, y).diff(x, y) - f(x, y).diff(y, x) == 0
    assert Derivative(f(x, y), x, y).args[1:] == (x, y)
    assert Derivative(f(x, y), y, x).args[1:] == (y, x)
    assert (Derivative(f(x, y), x, y) - Derivative(f(x, y), y, x)).doit() == 0


def test_suppressed_evaluation():
    a = sin(0, evaluate=False)
    assert a != 0
    assert a.func is sin
    assert a.args == (0,)


def test_function_evalf():
    def eq(a, b, eps):
        return abs(a - b) < eps
    assert eq(sin(1).evalf(15), Float("0.841470984807897"), 1e-13)
    assert eq(
        sin(2).evalf(25), Float("0.9092974268256816953960199", 25), 1e-23)
    assert eq(sin(1 + I).evalf(
        15), Float("1.29845758141598") + Float("0.634963914784736")*I, 1e-13)
    assert eq(exp(1 + I).evalf(15), Float(
        "1.46869393991588") + Float("2.28735528717884239")*I, 1e-13)
    assert eq(exp(-0.5 + 1.5*I).evalf(15), Float(
        "0.0429042815937374") + Float("0.605011292285002")*I, 1e-13)
    assert eq(log(pi + sqrt(2)*I).evalf(
        15), Float("1.23699044022052") + Float("0.422985442737893")*I, 1e-13)
    assert eq(cos(100).evalf(15), Float("0.86231887228768"), 1e-13)


def test_extensibility_eval():
    class MyFunc(Function):
        @classmethod
        def eval(cls, *args):
            return (0, 0, 0)
    assert MyFunc(0) == (0, 0, 0)


def test_function_non_commutative():
    x = Symbol('x', commutative=False)
    assert f(x).is_commutative is False
    assert sin(x).is_commutative is False
    assert exp(x).is_commutative is False
    assert log(x).is_commutative is False
    assert f(x).is_complex is False
    assert sin(x).is_complex is False
    assert exp(x).is_complex is False
    assert log(x).is_complex is False


def test_function_complex():
    x = Symbol('x', complex=True)
    assert f(x).is_commutative is True
    assert sin(x).is_commutative is True
    assert exp(x).is_commutative is True
    assert log(x).is_commutative is True
    assert f(x).is_complex is True
    assert sin(x).is_complex is True
    assert exp(x).is_complex is True
    assert log(x).is_complex is True


def test_function__eval_nseries():
    n = Symbol('n')

    assert sin(x)._eval_nseries(x, 2, None) == x + O(x**2)
    assert sin(x + 1)._eval_nseries(x, 2, None) == x*cos(1) + sin(1) + O(x**2)
    assert sin(pi*(1 - x))._eval_nseries(x, 2, None) == pi*x + O(x**2)
    assert acos(1 - x**2)._eval_nseries(x, 2, None) == sqrt(2)*x + O(x**2)
    assert polygamma(n, x + 1)._eval_nseries(x, 2, None) == \
        polygamma(n, 1) + polygamma(n + 1, 1)*x + O(x**2)
    raises(PoleError, lambda: sin(1/x)._eval_nseries(x, 2, None))
    raises(PoleError, lambda: acos(1 - x)._eval_nseries(x, 2, None))
    raises(PoleError, lambda: acos(1 + x)._eval_nseries(x, 2, None))
    assert loggamma(1/x)._eval_nseries(x, 0, None) == \
        log(x)/2 - log(x)/x - 1/x + O(1, x)
    assert loggamma(log(1/x)).nseries(x, n=1, logx=y) == loggamma(-y)

    # issue 6725:
    assert expint(S(3)/2, -x)._eval_nseries(x, 5, None) == \
        2 - 2*sqrt(pi)*sqrt(-x) - 2*x - x**2/3 - x**3/15 - x**4/84 + O(x**5)
    assert sin(sqrt(x))._eval_nseries(x, 3, None) == \
        sqrt(x) - x**(S(3)/2)/6 + x**(S(5)/2)/120 + O(x**3)


def test_doit():
    n = Symbol('n', integer=True)
    f = Sum(2 * n * x, (n, 1, 3))
    d = Derivative(f, x)
    assert d.doit() == 12
    assert d.doit(deep=False) == Sum(2*n, (n, 1, 3))


def test_evalf_default():
    from sympy.functions.special.gamma_functions import polygamma
    assert type(sin(4.0)) == Float
    assert type(re(sin(I + 1.0))) == Float
    assert type(im(sin(I + 1.0))) == Float
    assert type(sin(4)) == sin
    assert type(polygamma(2.0, 4.0)) == Float
    assert type(sin(Rational(1, 4))) == sin


def test_issue_5399():
    args = [x, y, S(2), S.Half]

    def ok(a):
        """Return True if the input args for diff are ok"""
        if not a:
            return False
        if a[0].is_Symbol is False:
            return False
        s_at = [i for i in range(len(a)) if a[i].is_Symbol]
        n_at = [i for i in range(len(a)) if not a[i].is_Symbol]
        # every symbol is followed by symbol or int
        # every number is followed by a symbol
        return (all(a[i + 1].is_Symbol or a[i + 1].is_Integer
            for i in s_at if i + 1 < len(a)) and
            all(a[i + 1].is_Symbol
            for i in n_at if i + 1 < len(a)))
    eq = x**10*y**8
    for a in subsets(args):
        for v in variations(a, len(a)):
            if ok(v):
                noraise = eq.diff(*v)
            else:
                raises(ValueError, lambda: eq.diff(*v))


def test_derivative_numerically():
    from random import random
    z0 = random() + I*random()
    assert abs(Derivative(sin(x), x).doit_numerically(z0) - cos(z0)) < 1e-15


def test_fdiff_argument_index_error():
    from sympy.core.function import ArgumentIndexError

    class myfunc(Function):
        nargs = 1  # define since there is no eval routine

        def fdiff(self, idx):
            raise ArgumentIndexError
    mf = myfunc(x)
    assert mf.diff(x) == Derivative(mf, x)
    raises(TypeError, lambda: myfunc(x, x))


def test_deriv_wrt_function():
    x = f(t)
    xd = diff(x, t)
    xdd = diff(xd, t)
    y = g(t)
    yd = diff(y, t)

    assert diff(x, t) == xd
    assert diff(2 * x + 4, t) == 2 * xd
    assert diff(2 * x + 4 + y, t) == 2 * xd + yd
    assert diff(2 * x + 4 + y * x, t) == 2 * xd + x * yd + xd * y
    assert diff(2 * x + 4 + y * x, x) == 2 + y
    assert (diff(4 * x**2 + 3 * x + x * y, t) == 3 * xd + x * yd + xd * y +
            8 * x * xd)
    assert (diff(4 * x**2 + 3 * xd + x * y, t) == 3 * xdd + x * yd + xd * y +
            8 * x * xd)
    assert diff(4 * x**2 + 3 * xd + x * y, xd) == 3
    assert diff(4 * x**2 + 3 * xd + x * y, xdd) == 0
    assert diff(sin(x), t) == xd * cos(x)
    assert diff(exp(x), t) == xd * exp(x)
    assert diff(sqrt(x), t) == xd / (2 * sqrt(x))


def test_diff_wrt_value():
    assert Expr()._diff_wrt is False
    assert x._diff_wrt is True
    assert f(x)._diff_wrt is True
    assert Derivative(f(x), x)._diff_wrt is True
    assert Derivative(x**2, x)._diff_wrt is False


def test_diff_wrt():
    fx = f(x)
    dfx = diff(f(x), x)
    ddfx = diff(f(x), x, x)

    assert diff(sin(fx) + fx**2, fx) == cos(fx) + 2*fx
    assert diff(sin(dfx) + dfx**2, dfx) == cos(dfx) + 2*dfx
    assert diff(sin(ddfx) + ddfx**2, ddfx) == cos(ddfx) + 2*ddfx
    assert diff(fx**2, dfx) == 0
    assert diff(fx**2, ddfx) == 0
    assert diff(dfx**2, fx) == 0
    assert diff(dfx**2, ddfx) == 0
    assert diff(ddfx**2, dfx) == 0

    assert diff(fx*dfx*ddfx, fx) == dfx*ddfx
    assert diff(fx*dfx*ddfx, dfx) == fx*ddfx
    assert diff(fx*dfx*ddfx, ddfx) == fx*dfx

    assert diff(f(x), x).diff(f(x)) == 0
    assert (sin(f(x)) - cos(diff(f(x), x))).diff(f(x)) == cos(f(x))

    assert diff(sin(fx), fx, x) == diff(sin(fx), x, fx)

    # Chain rule cases
    assert f(g(x)).diff(x) == \
        Subs(Derivative(f(x), x), (x,), (g(x),))*Derivative(g(x), x)
    assert diff(f(g(x), h(x)), x) == \
        Subs(Derivative(f(y, h(x)), y), (y,), (g(x),))*Derivative(g(x), x) + \
        Subs(Derivative(f(g(x), y), y), (y,), (h(x),))*Derivative(h(x), x)
    assert f(
        sin(x)).diff(x) == Subs(Derivative(f(x), x), (x,), (sin(x),))*cos(x)

    assert diff(f(g(x)), g(x)) == Subs(Derivative(f(x), x), (x,), (g(x),))


def test_diff_wrt_func_subs():
    assert f(g(x)).diff(x).subs(g, Lambda(x, 2*x)).doit() == f(2*x).diff(x)


def test_diff_wrt_not_allowed():
    raises(ValueError, lambda: diff(sin(x**2), x**2))
    raises(ValueError, lambda: diff(exp(x*y), x*y))
    raises(ValueError, lambda: diff(1 + x, 1 + x))


def test_klein_gordon_lagrangian():
    m = Symbol('m')
    phi = f(x, t)

    L = -(diff(phi, t)**2 - diff(phi, x)**2 - m**2*phi**2)/2
    eqna = Eq(
        diff(L, phi) - diff(L, diff(phi, x), x) - diff(L, diff(phi, t), t), 0)
    eqnb = Eq(diff(phi, t, t) - diff(phi, x, x) + m**2*phi, 0)
    assert eqna == eqnb


def test_sho_lagrangian():
    m = Symbol('m')
    k = Symbol('k')
    x = f(t)

    L = m*diff(x, t)**2/2 - k*x**2/2
    eqna = Eq(diff(L, x), diff(L, diff(x, t), t))
    eqnb = Eq(-k*x, m*diff(x, t, t))
    assert eqna == eqnb

    assert diff(L, x, t) == diff(L, t, x)
    assert diff(L, diff(x, t), t) == m*diff(x, t, 2)
    assert diff(L, t, diff(x, t)) == -k*x + m*diff(x, t, 2)


def test_straight_line():
    F = f(x)
    Fd = F.diff(x)
    L = sqrt(1 + Fd**2)
    assert diff(L, F) == 0
    assert diff(L, Fd) == Fd/sqrt(1 + Fd**2)


def test_sort_variable():
    vsort = Derivative._sort_variables

    assert vsort((x, y, z)) == [x, y, z]
    assert vsort((h(x), g(x), f(x))) == [f(x), g(x), h(x)]
    assert vsort((z, y, x, h(x), g(x), f(x))) == [x, y, z, f(x), g(x), h(x)]
    assert vsort((x, f(x), y, f(y))) == [x, f(x), y, f(y)]
    assert vsort((y, x, g(x), f(x), z, h(x), y, x)) == \
        [x, y, f(x), g(x), z, h(x), x, y]
    assert vsort((z, y, f(x), x, f(x), g(x))) == [y, z, f(x), x, f(x), g(x)]
    assert vsort((z, y, f(x), x, f(x), g(x), z, z, y, x)) == \
        [y, z, f(x), x, f(x), g(x), x, y, z, z]


def test_unhandled():
    class MyExpr(Expr):
        def _eval_derivative(self, s):
            if not s.name.startswith('xi'):
                return self
            else:
                return None

    expr = MyExpr(x, y, z)
    assert diff(expr, x, y, f(x), z) == Derivative(expr, f(x), z)
    assert diff(expr, f(x), x) == Derivative(expr, f(x), x)


def test_issue_4711():
    x = Symbol("x")
    assert Symbol('f')(x) == f(x)


def test_nfloat():
    from sympy.core.basic import _aresame
    from sympy.polys.rootoftools import rootof

    x = Symbol("x")
    eq = x**(S(4)/3) + 4*x**(S(1)/3)/3
    assert _aresame(nfloat(eq), x**(S(4)/3) + (4.0/3)*x**(S(1)/3))
    assert _aresame(nfloat(eq, exponent=True), x**(4.0/3) + (4.0/3)*x**(1.0/3))
    eq = x**(S(4)/3) + 4*x**(x/3)/3
    assert _aresame(nfloat(eq), x**(S(4)/3) + (4.0/3)*x**(x/3))
    big = 12345678901234567890
    # specify precision to match value used in nfloat
    Float_big = Float(big, 15)
    assert _aresame(nfloat(big), Float_big)
    assert _aresame(nfloat(big*x), Float_big*x)
    assert _aresame(nfloat(x**big, exponent=True), x**Float_big)
    assert nfloat({x: sqrt(2)}) == {x: nfloat(sqrt(2))}
    assert nfloat({sqrt(2): x}) == {sqrt(2): x}
    assert nfloat(cos(x + sqrt(2))) == cos(x + nfloat(sqrt(2)))

    # issue 6342
    f = S('x*lamda + lamda**3*(x/2 + 1/2) + lamda**2 + 1/4')
    assert not any(a.free_symbols for a in solve(f.subs(x, -0.139)))

    # issue 6632
    assert nfloat(-100000*sqrt(2500000001) + 5000000001) == \
        9.99999999800000e-11

    # issue 7122
    eq = cos(3*x**4 + y)*rootof(x**5 + 3*x**3 + 1, 0)
    assert str(nfloat(eq, exponent=False, n=1)) == '-0.7*cos(3.0*x**4 + y)'


def test_issue_7068():
    from sympy.abc import a, b, f
    y1 = Dummy('y')
    y2 = Dummy('y')
    func1 = f(a + y1 * b)
    func2 = f(a + y2 * b)
    func1_y = func1.diff(y1)
    func2_y = func2.diff(y2)
    assert func1_y != func2_y
    z1 = Subs(f(a), a, y1)
    z2 = Subs(f(a), a, y2)
    assert z1 != z2


def test_issue_7231():
    from sympy.abc import a
    ans1 = f(x).series(x, a)
    _xi_1 = ans1.atoms(Dummy).pop()
    res = (f(a) + (-a + x)*Subs(Derivative(f(_xi_1), _xi_1), (_xi_1,), (a,)) +
           (-a + x)**2*Subs(Derivative(f(_xi_1), _xi_1, _xi_1), (_xi_1,), (a,))/2 +
           (-a + x)**3*Subs(Derivative(f(_xi_1), _xi_1, _xi_1, _xi_1),
                            (_xi_1,), (a,))/6 +
           (-a + x)**4*Subs(Derivative(f(_xi_1), _xi_1, _xi_1, _xi_1, _xi_1),
                            (_xi_1,), (a,))/24 +
           (-a + x)**5*Subs(Derivative(f(_xi_1), _xi_1, _xi_1, _xi_1, _xi_1, _xi_1),
                            (_xi_1,), (a,))/120 + O((-a + x)**6, (x, a)))
    assert res == ans1
    ans2 = f(x).series(x, a)
    assert res == ans2


def test_issue_7687():
    from sympy.core.function import Function
    from sympy.abc import x
    f = Function('f')(x)
    ff = Function('f')(x)
    match_with_cache = ff.matches(f)
    assert isinstance(f, type(ff))
    clear_cache()
    ff = Function('f')(x)
    assert isinstance(f, type(ff))
    assert match_with_cache == ff.matches(f)


def test_issue_7688():
    from sympy.core.function import Function, UndefinedFunction

    f = Function('f')  # actually an UndefinedFunction
    clear_cache()
    class A(UndefinedFunction):
        pass
    a = A('f')
    assert isinstance(a, type(f))


def test_mexpand():
    from sympy.abc import x
    assert _mexpand(None) is None
    assert _mexpand(1) is S.One
    assert _mexpand(x*(x + 1)**2) == (x*(x + 1)**2).expand()
