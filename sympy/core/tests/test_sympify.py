from sympy import Symbol, exp, Integer, Real, sin, cos, log, Poly, Lambda, \
        Function, I, S, sqrt,  raises, srepr
from sympy.abc import x, y
from sympy.core.sympify import sympify, _sympify, SympifyError
from sympy.core.decorators import _sympifyit

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
    assert sympify("x^3") == x**3
    assert sympify("1/2") == Integer(1)/2

    raises(SympifyError, "_sympify('x**3')")
    raises(SympifyError, "_sympify('1/2')")

def test_sympify_bool():
    """Test that sympify accepts boolean values
    and that output leaves them unchanged"""
    assert sympify(True) == True
    assert sympify(False)== False

def test_sympify4():
    class A:
        def _sympy_(self):
            return Symbol("x")

    a = A()

    assert _sympify(a)**3== x**3
    assert sympify(a)**3 == x**3
    assert a == x

def test_sympify_text():
    assert sympify('some') == Symbol('some')
    assert sympify('core') == Symbol('core')

    assert sympify('True') == True
    assert sympify('False') == False

    assert sympify('Poly') == Poly
    assert sympify('sin') == sin

def test_sympify_function():
    assert sympify('factor(x**2-1, x)') == -(1-x)*(x+1)
    assert sympify('sin(pi/2)*cos(pi)') == -Integer(1)

def test_sympify_poly():
    p = Poly(x**2+x+1, x)

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
    assert sympify('lambda : 1') == Lambda(x, 1)
    assert sympify('lambda x: 2*x') == Lambda(x, 2*x)
    assert sympify('lambda x, y: 2*x+y') == Lambda([x, y], 2*x+y)

    raises(SympifyError, "_sympify('lambda : 1')")

def test_sympify_raises():
    raises(SympifyError, 'sympify("fx)")')


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
    raises(SympifyError, "_sympify('1')")
    raises(SympifyError, "_sympify([1,2,3])")


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

    raises(SympifyError, "add_raises(x, '1')")

def test_int_float():
    class F1_1(object):
        def __float__(self):
            return 1.1

    class F1_1b(object):
        """
        This class is still a float, even though it also implements __int__().
        """
        def __float__(self):
            return 1.1

        def __int__(self):
            return 1

    class F1_1c(object):
        """
        This class is still a float, because it implements _sympy_()
        """
        def __float__(self):
            return 1.1

        def __int__(self):
            return 1

        def _sympy_(self):
            return Real(1.1)

    class I5(object):
        def __int__(self):
            return 5

    class I5b(object):
        """
        This class implements both __int__() and __float__(), so it will be
        treated as Real in SymPy. One could change this behavior, by using
        float(a) == int(a), but deciding that integer-valued floats represent
        exact numbers is arbitrary and often not correct, so we do not do it.
        If, in the future, we decide to do it anyway, the tests for I5b need to
        be changed.
        """
        def __float__(self):
            return 5.0

        def __int__(self):
            return 5

    class I5c(object):
        """
        This class implements both __int__() and __float__(), but also
        a _sympy_() method, so it will be Integer.
        """
        def __float__(self):
            return 5.0

        def __int__(self):
            return 5

        def _sympy_(self):
            return Integer(5)

    i5 = I5()
    i5b = I5b()
    i5c = I5c()
    f1_1 = F1_1()
    f1_1b = F1_1b()
    f1_1c = F1_1c()
    assert sympify(i5) == 5
    assert isinstance(sympify(i5), Integer)
    assert sympify(i5b) == 5
    assert isinstance(sympify(i5b), Real)
    assert sympify(i5c) == 5
    assert isinstance(sympify(i5c), Integer)
    assert abs(sympify(f1_1) - 1.1) < 1e-5
    assert abs(sympify(f1_1b) - 1.1) < 1e-5
    assert abs(sympify(f1_1c) - 1.1) < 1e-5

    assert _sympify(i5) == 5
    assert isinstance(_sympify(i5), Integer)
    assert _sympify(i5b) == 5
    assert isinstance(_sympify(i5b), Real)
    assert _sympify(i5c) == 5
    assert isinstance(_sympify(i5c), Integer)
    assert abs(_sympify(f1_1) - 1.1) < 1e-5
    assert abs(_sympify(f1_1b) - 1.1) < 1e-5
    assert abs(_sympify(f1_1c) - 1.1) < 1e-5


def test_issue1034():
    a = sympify('Integer(4)')

    assert a == Integer(4)
    assert a.is_Integer

def test_issue883():
    a = [3,2.0]
    assert sympify(a) == [Integer(3), Real(2.0)]
    assert sympify(tuple(a)) == (Integer(3), Real(2.0))
    assert sympify(set(a)) == set([Integer(3), Real(2.0)])

def test_S_sympify():
    assert S(1)/2 == sympify(1)/2
    assert (-2)**(S(1)/2) == sqrt(2)*I

def test_issue1689():
    assert srepr(S(1.0+0J)) == srepr(S(1.0)) == srepr(Real(1.0))
    assert srepr(Real(1)) != srepr(Real(1.0))

def test_issue1699_None():
    assert S(None) == None
