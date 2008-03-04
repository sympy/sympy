from sympy import Symbol, exp, Integer, Real, sin, cos, log, Polynomial, Lambda, Function, I
from sympy.abc import x, y
from sympy.core.sympify import sympify, _sympify, _sympifyit, SympifyError
import py
#from sympy.utilities.pytest import XFAIL


def test_439():
    v = sympify("exp(x)")
    x = Symbol("x")
    assert v == exp(x)
    assert type(v) == type(exp(x))
    assert str(type(v)) == str(type(exp(x)))

def test_sympify1():
    assert sympify("x") == Symbol("x")
    assert sympify("   x") == Symbol("x")
    assert sympify("   x   ") == Symbol("x")

def test_sympify2():
    class A:
        def _sympy_(self):
            return Symbol("x")**3

    a = A()

    assert _sympify(a)== x**3
    assert sympify(a) == x**3
    assert a == x**3

def test_sympify3():
    assert sympify("x**3") == x**3
    assert sympify("1/2") == Integer(1)/2

    py.test.raises(SympifyError, "_sympify('x**3')")
    py.test.raises(SympifyError, "_sympify('1/2')")

def test_sympify4():
    class A:
        def _sympy_(self):
            return Symbol("x")

    a = A()

    assert _sympify(a)**3== x**3
    assert sympify(a)**3 == x**3
    assert a == x


def test_sympify_poly():
    p = Polynomial(1+x+x**2)

    assert _sympify(p) is p
    assert sympify(p) is p

def test_sage():
    # how to effectivelly test for the _sage_() method without having SAGE
    # installed?
    assert hasattr(x, "_sage_")
    assert hasattr(Integer(3), "_sage_")
    assert hasattr(sin(x), "_sage_")
    assert hasattr(cos(x), "_sage_")
    assert hasattr(x**2, "_sage_")
    assert hasattr(x+y, "_sage_")
    assert hasattr(exp(x), "_sage_")
    assert hasattr(log(x), "_sage_")

def test_bug496():
    a_ = sympify("a_")
    _a = sympify("_a")

def test_lambda():
    x = Symbol('x')
    assert sympify('lambda : 1')==Lambda(x, 1)
    assert sympify('lambda x: 2*x')==Lambda(x, 2*x)

    py.test.raises(SympifyError, "_sympify('lambda : 1')")

def test_sympify_raises():
    py.test.raises(SympifyError, sympify, "fx)")


def test__sympify():
    x = Symbol('x')
    f = Function('f')

    # positive _sympify
    assert _sympify(x)      is x
    assert _sympify(f)      is f
    assert _sympify(1)      == Integer(1)
    assert _sympify(0.5)    == Real("0.5")
    assert _sympify(1+1j)   == 1 + I

    class A:
        def _sympy_(self):
            return Integer(5)

    a = A()
    assert _sympify(a)      == Integer(5)

    # negative _sympify
    py.test.raises(SympifyError, "_sympify('1')")
    py.test.raises(SympifyError, "_sympify([1,2,3])")


def test_sympifyit():
    x = Symbol('x')
    y = Symbol('y')

    @_sympifyit('b', NotImplemented)
    def add(a, b):
        return a+b

    assert add(x, 1)    == x+1
    assert add(x, 0.5)  == x+Real('0.5')
    assert add(x, y)    == x+y

    assert add(x, '1')  == NotImplemented


    @_sympifyit('b')
    def add_raises(a, b):
        return a+b

    assert add_raises(x, 1)     == x+1
    assert add_raises(x, 0.5)   == x+Real('0.5')
    assert add_raises(x, y)     == x+y

    py.test.raises(SympifyError, "add_raises(x, '1')")
