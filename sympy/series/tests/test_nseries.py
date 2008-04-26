from sympy import Symbol, Rational, ln, exp, log, sqrt, E, O
from sympy.abc import x, y, z
from sympy.utilities.pytest import XFAIL

def test_simple_1():
    assert x.nseries(x, 0, 5) == x
    assert y.nseries(x, 0, 5) == y
    assert (1/(x*y)).nseries(y, 0, 5) == 1/(x*y)
    assert Rational(3,4).nseries(x, 0, 5) == Rational(3,4)

def test_mul_0():
    assert (x*ln(x)).nseries(x, 0, 5) == x*ln(x)

def test_mul_1():
    assert (x*ln(2+x)).nseries(x, 0, 4) == x*log(2)+x**2/2-x**3/8+x**4/24+ \
            O(x**5)
    assert (x*ln(1+x)).nseries(x, 0, 4) == x**2- x**3/2 + x**4/3 + O(x**5)

def test_pow_0():
    assert (x**2).nseries(x, 0, 5) == x**2
    assert (1/x).nseries(x, 0, 5) == 1/x
    assert (1/x**2).nseries(x, 0, 5) == 1/x**2
    assert (x**(Rational(2,3))).nseries(x, 0, 5) == (x**(Rational(2,3)))
    assert (x**(Rational(3,2))).nseries(x, 0, 5) == (x**(Rational(3,2)))

def test_pow_1():
    assert ((1+x)**2).nseries(x, 0, 5) == 1+2*x+x**2

def test_geometric_1():
    assert (1/(1-x)).nseries(x, 0, 5) == 1+x+x**2+x**3+x**4+O(x**5)
    assert (x/(1-x)).nseries(x, 0, 5) == x+x**2+x**3+x**4+x**5+O(x**6)
    assert (x**3/(1-x)).nseries(x, 0, 5) == x**3+x**4+x**5+x**6+x**7+O(x**8)

def test_sqrt_1():
    assert sqrt(1+x).nseries(x, 0, 5) == 1+x/2-x**2/8+x**3/16-5*x**4/128+O(x**5)

def test_exp_1():
    assert exp(x).nseries(x, 0, 5) == 1+x+x**2/2+x**3/6+x**4/24 + O(x**5)
    assert exp(x).nseries(x, 0, 12) == 1+x+x**2/2+x**3/6+x**4/24+x**5/120+  \
            x**6/720+x**7/5040+x**8/40320+x**9/362880+x**10/3628800+  \
            x**11/39916800 + O(x**12)
    assert exp(1/x).nseries(x, 0, 5) == exp(1/x)
    assert exp(1/(1+x)).nseries(x, 0, 4) ==  \
            (E*(1-x-13*x**3/6+3*x**2/2)).expand() + O(x**4)
    assert exp(2+x).nseries(x, 0, 5) ==  \
            (exp(2)*(1+x+x**2/2+x**3/6+x**4/24)).expand() + O(x**5)

def test_exp_sqrt_1():
    assert exp(1+sqrt(x)).nseries(x, 0, 2) ==  \
        (exp(1)*(1+sqrt(x)+x/2+sqrt(x)*x/6)).expand() + O(x**2)

def test_power_x_x1():
    assert (exp(x*ln(x))).nseries(x, 0, 4) == \
            1+x*log(x)+x**2*log(x)**2/2+x**3*log(x)**3/6 + O(x**4*log(x)**4)

@XFAIL
def test_power_x_x2():
    assert (x**x).nseries(x, 0, 4) == \
            1+x*log(x)+x**2*log(x)**2/2+x**3*log(x)**3/6 + O(x**4*log(x)**4)
