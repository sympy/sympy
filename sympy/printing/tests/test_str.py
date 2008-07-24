from sympy import Symbol, Wild, Rational, Derivative, I, log, sqrt, exp, sin, abs,\
                  factorial, Lambda, O, cos, Function, WildFunction, sum, zeta
from sympy.polynomials import Polynomial
from sympy.polys.polynomial import Poly
from sympy.polys.rootfinding import RootsOf, RootOf, RootSum
from sympy.utilities.pytest import XFAIL

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
w = Symbol('w')

def test_Symbol():
    assert str(x)=="x"
    e = x
    assert str(e)=="x"

def test_Polynomial():
    f = Polynomial(x+2)
    assert str(f) == "2 + x"

def test_poly():
    assert str(Poly(0, x)) == "Poly(0, x)"

    assert str(Poly(1, x)) == "Poly(1, x)"
    assert str(Poly(x, x)) == "Poly(x, x)"

    assert str(Poly(2*x + 1, x)) == "Poly(2*x + 1, x)"
    assert str(Poly(2*x - 1, x)) == "Poly(2*x - 1, x)"

    assert str(Poly(-1, x)) == "Poly(-1, x)"
    assert str(Poly(-x, x)) == "Poly(-x, x)"

    assert str(Poly(-2*x + 1, x)) == "Poly(-2*x + 1, x)"
    assert str(Poly(-2*x - 1, x)) == "Poly(-2*x - 1, x)"

    assert str(Poly(x - 1, x, order='lex')) == "Poly(x - 1, x)"

    assert str(Poly(x**2 + 1 + y, x)) == "Poly(x**2 + 1 + y, x)"
    assert str(Poly(x**2 - 1 + y, x)) in [
            "Poly(x**2 - 1 + y, x)",
            "Poly(x**2 + (-1) + y, x)",
            ]

    assert str(Poly(-x*y*z + x*y - 1, x, y, z)) == "Poly(-x*y*z + x*y - 1, x, y, z)"

    assert str(Poly(-w*x**21*y**7*z + (1 + w)*z**3 - 2*x*z + 1, x, y, z)) == \
        "Poly(-w*x**21*y**7*z + (1 + w)*z**3 - 2*x*z + 1, x, y, z)"

    assert str(Poly(x*y*z**2 - 27*x, x, y, z, order='lex')) == \
        "Poly(x*y*z**2 - 27*x, x, y, z, order='lex')"
    assert str(Poly(x*y*z**2 - 27*x, x, y, z, order='grlex')) == \
        "Poly(x*y*z**2 - 27*x, x, y, z)"
    assert str(Poly(x*y*z**2 - 27*x, x, y, z, order='grevlex')) == \
        "Poly(x*y*z**2 - 27*x, x, y, z, order='grevlex')"

def test_RootsOf():
    f = Poly(x**17 + x - 1, x)

    assert str(RootsOf(f)) == "RootsOf(x**17 + x - 1, x)"

    assert str(RootOf(f, 0)) == "RootOf(x**17 + x - 1, x, index=0)"

    assert str(RootSum(Lambda(z, z**2), f)) == \
        "RootSum(Lambda(_z, _z**2), x**17 + x - 1, x)"

def test_order():
    assert str(O(x)) == "O(x)"
    assert str(O(x**2)) == "O(x**2)"

    assert str(O(x*y)) == "O(x*y, x, y)"

def test_pow():
    assert str(1/x) == "1/x"

def test_ordering():
    from sympy.core.methods import ArithMeths
    class CustomClass1(Basic, ArithMeths): pass
    class CustomClass2(Basic, ArithMeths): pass
    cc1 = CustomClass1(commutative=True)
    cc2 = CustomClass2(commutative=True)

    assert str(Rational(2)*cc1) == '2*CustomClass1()'
    assert str(cc1*Rational(2)) == '2*CustomClass1()'
    assert str(cc1*Real("1.5")) == '1.5*CustomClass1()'
    assert str(cc2*Rational(2)) == '2*CustomClass2()'
    assert str(cc2*Rational(2)*cc1) == '2*CustomClass1()*CustomClass2()'
    assert str(cc1*Rational(2)*cc2) == '2*CustomClass1()*CustomClass2()'

def test_abs():
    assert str(abs(x)) == "abs(x)"
    assert str(abs(Rational(1,6))) == "1/6"
    assert str(abs(Rational(-1,6))) == "1/6"

def test_concrete():
    assert str(sum(cos(3*z), (z, x, y))) == "Sum(cos(3*z), (z, x, y))"

def test_zeta():
    assert str(zeta(3)) == "zeta(3)"

def test_factorials():
    n = Symbol('n', integer=True)
    assert str(factorial(-2)) == "0"
    assert str(factorial(0)) == "1"
    assert str(factorial(7)) == "5040"
    assert str(factorial(n)) == "n!"
    assert str(factorial(2*n)) == "(2*n)!"

def test_add_eval():
    e = x+y+x+y
    s1 = str(e)
    s2 = str(e)
    assert s1 == s2

def test_function_str():
    f = Function('f')
    fx= f(x)
    w = WildFunction('w')
    assert str(fx)  == "f(x)"
    assert str(w)   == "w_"

@XFAIL
def test_unapplied_function_str():
    f = Function('f')
    assert repr(f)     == "Function('f')"   # this does not work
    assert str(f)      == "f"               # this does not work

def test_Lambda():
    d = Symbol('d', dummy=True)
    assert str(Lambda(d, d**2)) == "Lambda(_d, _d**2)"

def test_suppressed_evaluation():
    a = sin(0, evaluate=False)
    assert str(a) == "sin(0)"

def test_Rational():
    n1 = Rational(1,4)
    n2 = Rational(1,3)
    n3 = Rational(2,4)
    n4 = Rational(2,-4)
    n5 = Rational(0)
    n6 = Rational(1)
    n7 = Rational(3)
    n8 = Rational(-3)
    assert str(n1*n2) == "1/12"
    assert str(n1*n2) == "1/12"
    assert str(n3) == "1/2"
    assert str(n1*n3) == "1/8"
    assert str(n1+n3) == "3/4"
    assert str(n1+n2) == "7/12"
    assert str(n1+n4) == "-1/4"
    assert str(n4*n4) == "1/4"
    assert str(n4+n2) == "-1/6"
    assert str(n4+n5) == "-1/2"
    assert str(n4*n5) == "0"
    assert str(n3+n4) == "0"
    assert str(n1**n7) == "1/64"
    assert str(n2**n7) == "1/27"
    assert str(n2**n8) == "27"
    assert str(n7**n8) == "1/27"
    assert str(Rational("-25")) == "-25"
    assert str(Rational("25/7")) == "25/7"
    assert str(Rational("-123/569")) == "-123/569"
    assert str(Rational("1.25")) == "5/4"
    assert str(Rational("-2.6e-2")) == "-13/500"
    assert str(Rational("0.1[23]")) == "61/495"
    assert str(Rational("5.1[666]")) == "31/6"
    assert str(Rational("-5.1[666]")) == "-31/6"
    assert str(Rational("0.[9]")) == "1"
    assert str(Rational("-0.[9]")) == "-1"

    assert str(Rational(1,4) ** Rational(1,2)) == "1/2"
    assert str(Rational(1,36) ** Rational(1,2)) == "1/6"

    assert str((123**25) ** Rational(1,25)) == "123"
    assert str((123**25+1)**Rational(1,25)) != "123"
    assert str((123**25-1)**Rational(1,25)) != "123"
    assert str((123**25-1)**Rational(1,25)) != "122"

    assert str(Rational(81,36)**(Rational(3,2))) == "27/8"
    assert str(Rational(81,36)**(-Rational(3,2))) == "8/27"

    assert str((-4)**Rational(1,2)) == str(2*I)
    assert str(2**Rational(1,10**10)) == "2**(1/10000000000)"

def test_poly_str():
    #if any of these tests fails, it can still be correct, just the terms can
    #be in a different order. That happens for example when we change the
    #hash algorithm. If it is correct, just add another item in the list [] of
    #correct results.
    assert str((2*x-(7*x**2 - 2) + 3*y)) in ["2 - 7*x**2 + 2*x + 3*y",
            "2 + 3*y + 2*x - 7*x**2", "2 + 3*y - 7*x**2 + 2*x",
            "3*y + 2*x + 2 - 7*x**2", "2 + 2*x + 3*y - 7*x**2"]
    assert str(x-y) in ["x - y", "-y + x"]
    assert str(2+-x) in ["2 - x", "-x + 2"]
    assert str(x-2) in ["x - 2", "(-2) + x", "-2 + x"]
    assert str(x-y-z-w) in ["x - y - z - w","-w - y - z + x","x - w - y - z",
                            "-w + x - y - z","-z - w - y + x","-y + x - w - z",
                            "-y - z - w + x"]
    assert str(x-y-z-w) in [
            "-w - y - z + x","x - w - y - z","-w + x - z - y","-y - w - z + x",
            "-y + x - z - w","-y + x - w - z","-w + x - y - z","-z - w - y +x",
            "-y - z - w + x"]
    assert str(x-z*y**2*z*w) in ["-z**2*y**2*w + x", "x - w*y**2*z**2",
            "-y**2*z**2*w + x","x - w*z**2*y**2","x - y**2*z**2*w",
            "x - y**2*w*z**2","x - z**2*y**2*w","-w*z**2*y**2 + x",
            "-w*y**2*z**2 + x","x - z**2*w*y**2"]

def test_bug1():
    assert str(x-1*y*x*y) in ["x - x*y**2", "-x*y**2 + x"]

def test_bug2():
    e = x-y
    a = str(e)
    b = str(e)
    assert a == b

def test_bug3():
    x = Symbol("x")
    e = sqrt(x)
    assert str(e) == "x**(1/2)"

def test_bug4():
    e = -2*sqrt(x)-y/sqrt(x)/2
    assert str(e) not in ["(-2)*x**1/2(-1/2)*x**(-1/2)*y",
            "-2*x**1/2(-1/2)*x**(-1/2)*y","-2*x**1/2-1/2*x**-1/2*w"]
    assert str(e) in ["-2*x**(1/2) - 1/2*x**(-1/2)*y", "-2*x**(1/2) - 1/2*y*x**(-1/2)",
                      "-1/2*x**(-1/2)*y - 2*x**(1/2)", "-1/2*y*x**(-1/2) - 2*x**(1/2)",
                      "-2*x**(1/2) - y/(2*x**(1/2))"]

def test_Derivative():
    e = Derivative(x**2, x, evaluate=False)
    assert str(e) == "D(x**2, x)"

    e = Derivative(x**2/y, x, y, evaluate=False)
    assert str(e) == "D(1/y*x**2, x, y)"

def test_x_div_y():
    assert str(x/y) == "x/y"
    assert str(y/x) == "y/x"

def test_ordering():
    assert str(sin(x).series(x, 0, 15)) == "x - 1/6*x**3 + 1/120*x**5 - 1/5040*x**7 + 1/362880*x**9 - 1/39916800*x**11 + 1/6227020800*x**13 + O(x**15)"

def test_wild_str():
    # Check expressions containing Wild not causing infinite recursion
    w = Wild('x')
    assert str(w + 1) == str(x + 1)
    assert str(exp(2**w) + 5) == str(exp(2**x) + 5)
    assert str(3*w + 1) == str(3*x + 1)
    assert str(1/w + 1) == str(1/x + 1)
    assert str(w**2 + 1) == str(x**2 + 1)
    assert str(1/(1-w)) == str(1/(1-x))
