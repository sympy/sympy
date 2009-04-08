from sympy.utilities.pytest import XFAIL
from sympy import (symbols, lambdify, sqrt, sin, cos, pi, atan, Rational, Real,
        Matrix, Lambda, exp, Integral, oo)
from sympy import mpmath
import math, sympy

# high precision output of sin(0.2*pi) is used to detect if precision is lost unwanted
mpmath.mp.dps = 50
sin02 = mpmath.mpf("0.19866933079506121545941262711838975037020672954020")

x,y,z = symbols('xyz')

#================== Test different arguments ==============
def test_no_args():
    f = lambdify([], 1)
    try:
        f(-1)
        assert False
    except TypeError:
        pass
    assert f() == 1

def test_single_arg():
    f = lambdify(x, 2*x)
    assert f(1) == 2

def test_list_args():
    f = lambdify([x,y], x+y)
    assert f(1,2) == 3

def test_str_args():
    f = lambdify('x,y,z', 'z,y,x')
    assert f(3,2,1) == (1,2,3)
    assert f(1.0,2.0,3.0) == (3.0,2.0,1.0)
    # make sure correct number of args required
    try:
        f(0)
        assert False
    except TypeError:
        pass

def test_own_namespace():
    myfunc = lambda x:1
    f = lambdify(x, sin(x), {"sin":myfunc})
    assert f(0.1) == 1
    assert f(100) == 1

def test_own_module():
    f = lambdify(x, sin(x), math)
    assert f(0)==0.0
    f = lambdify(x, sympy.ceiling(x), math)
    try:
        f(4.5)
        assert False
    except NameError:
        pass

def test_bad_args():
    try:
        # no vargs given
        f = lambdify(1)
        assert False
    except TypeError:
        pass
    try:
        # same with vector exprs
        f = lambdify([1,2])
        assert False
    except TypeError:
        pass

#================== Test different modules ================
def test_sympy_lambda():
    f = lambdify(x, sin(x), "sympy")
    assert f(x) is sin(x)
    prec = 1e-15
    assert -prec < f(Rational(1,5)).evalf() - Real(str(sin02)) < prec
    try:
        # arctan is in numpy module and should not be available
        f = lambdify(x, arctan(x), "sympy")
        assert False
    except NameError:
        pass

def test_math_lambda():
    f = lambdify(x, sin(x), "math")
    prec = 1e-15
    assert -prec < f(0.2) - sin02 < prec
    try:
        f(x) # if this succeeds, it can't be a python math function
        assert False
    except ValueError:
        pass

def test_mpmath_lambda():
    f = lambdify(x, sin(x), "mpmath")
    prec = 1e-49 # mpmath precision is around 50 decimal places
    assert -prec < f(mpmath.mpf("0.2")) - sin02 < prec
    try:
        f(x) # if this succeeds, it can't be a mpmath function
        assert False
    except TypeError:
        pass

@XFAIL
def test_number_precision():
    f = lambdify(x, sin02, "mpmath")
    prec = 1e-49 # mpmath precision is around 50 decimal places
    assert -prec < f(0) - sin02 < prec

#================== Test Translations =====================
# We can only check if all translated functions are valid. It has to be checked
# by hand if they are complete.

def test_math_transl():
    from sympy.utilities.lambdify import MATH_TRANSLATIONS
    for sym, mat in MATH_TRANSLATIONS.iteritems():
        assert sym in sympy.functions.__dict__
        assert mat in math.__dict__

def test_mpmath_transl():
    from sympy.utilities.lambdify import MPMATH_TRANSLATIONS
    for sym, mat in MPMATH_TRANSLATIONS.iteritems():
        assert sym in sympy.functions.__dict__ or sym == 'Matrix'
        assert mat in mpmath.__dict__

#================== Test some functions ===================
def test_exponentiation():
    f = lambdify(x, x**2)
    assert f(-1) == 1
    assert f(0) == 0
    assert f(1) == 1
    assert f(-2) == 4
    assert f(2) == 4
    assert f(2.5) == 6.25

def test_sqrt():
    f = lambdify(x, sqrt(x))
    assert f(0) == 0.0
    assert f(1) == 1.0
    assert f(4) == 2.0
    assert abs(f(2) - 1.414) < 0.001
    assert f(6.25) == 2.5
    try:
        f(-1)
        assert False
    except ValueError: pass

def test_trig():
    f = lambdify([x], [cos(x),sin(x)])
    d = f(pi)
    prec = 1e-11
    assert -prec < d[0]+1 < prec
    assert -prec < d[1] < prec
    d = f(3.14159)
    prec = 1e-5
    assert -prec < d[0]+1 < prec
    assert -prec < d[1] < prec

#================== Test vectors ==========================
def test_vector_simple():
    f = lambdify((x,y,z), (z,y,x))
    assert f(3,2,1) == (1,2,3)
    assert f(1.0,2.0,3.0) == (3.0,2.0,1.0)
    # make sure correct number of args required
    try:
        f(0)
        assert False
    except TypeError: pass

def test_vector_discontinuous():
    f = lambdify(x, (-1/x, 1/x))
    try:
        f(0)
        assert False
    except ZeroDivisionError: pass
    assert f(1) == (-1.0, 1.0)
    assert f(2) == (-0.5, 0.5)
    assert f(-2) == (0.5, -0.5)

def test_trig_symbolic():
    f = lambdify([x], [cos(x),sin(x)])
    d = f(pi)
    assert abs(d[0]+1) < 0.0001
    assert abs(d[1]-0) < 0.0001

def test_trig_float():
    f = lambdify([x], [cos(x),sin(x)])
    d = f(3.14159)
    assert abs(d[0]+1) < 0.0001
    assert abs(d[1]-0) < 0.0001

def test_docs():
    f = lambdify(x, x**2)
    assert f(2) == 4
    f = lambdify([x,y,z], [z,y,x])
    assert f(1, 2, 3) == [3, 2, 1]
    f = lambdify(x, sqrt(x))
    assert f(4) == 2.0
    f = lambdify((x,y), sin(x*y)**2)
    assert f(0, 5) == 0

def test_math():
    f = lambdify((x, y), sin(x), modules="math")
    assert f(0, 5) == 0

def test_sin():
    f = lambdify(x, sin(x)**2)
    assert isinstance(f(2), float)
    f = lambdify(x, sin(x)**2, modules="math")
    assert isinstance(f(2), float)

def test_matrix():
    A = Matrix([[x, x*y], [sin(z)+4, x**z]])
    sol = Matrix([[1, 2], [sin(3)+4, 1]])
    f = lambdify((x,y,z), A, modules="sympy")
    assert f(1,2,3) == sol
    f = lambdify((x,y,z), (A, [A]), modules="sympy")
    assert f(1,2,3) == (sol,[sol])

def test_integral():
    f = Lambda(x, exp(-x**2))
    l = lambdify(x, Integral(f(x), (x, -oo, oo)), modules="sympy")
    assert l(x) == Integral(exp(-x**2), (x, -oo, oo))
