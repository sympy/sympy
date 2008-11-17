from sympy import Symbol, Rational, Order, C, exp, ln, log, O, var, nan, pi, S
from sympy.utilities.pytest import XFAIL
from sympy.abc import w, x, y, z

def test_caching_bug():
    #needs to be a first test, so that all caches are clean
    #cache it
    e = O(w)
    #and test that this won't raise an exception
    f = O(w**(-1/x/log(3)*log(5)), w)


def test_simple_1():
    o = Rational(0)
    assert Order(2*x) == Order(x)
    assert Order(x)*3 == Order(x)
    assert -28*Order(x) == Order(x)
    assert Order(-23) == Order(1)
    assert Order(exp(x)) == Order(1,x)
    assert Order(exp(1/x)).expr == exp(1/x)
    assert Order(x*exp(1/x)).expr == x*exp(1/x)
    assert Order(x**(o/3)).expr == x**(o/3)
    assert Order(x**(5*o/3)).expr == x**(5*o/3)

def test_simple_2():
    assert Order(2*x)*x == Order(x**2)
    assert Order(2*x)/x == Order(1,x)
    assert Order(2*x)*x*exp(1/x) == Order(x**2*exp(1/x))
    assert (Order(2*x)*x*exp(1/x)/ln(x)**3).expr == x**2*exp(1/x)*ln(x)**-3

def test_simple_3():
    assert Order(x)+x == Order(x)
    assert Order(x)+2 == 2+Order(x)
    assert Order(x)+x**2 == Order(x)
    assert Order(x)+1/x == 1/x+Order(x)
    assert Order(1/x)+1/x**2 == 1/x**2+Order(1/x)
    assert Order(x)+exp(1/x) == Order(x)+exp(1/x)

def test_simple_4():
    assert Order(x)**2 == Order(x**2)
    assert Order(x**3)**-2 == Order(x**-6)

def test_simple_5():
    assert Order(x)+Order(x**2) == Order(x)
    assert Order(x)+Order(x**-2) == Order(x**-2)
    assert Order(x)+Order(1/x) == Order(1/x)

def test_simple_6():
    assert Order(x)-Order(x) == Order(x)
    assert Order(x)+Order(1) == Order(1)
    assert Order(x)+Order(x**2) == Order(x)
    assert Order(1/x)+Order(1) == Order(1/x)
    assert Order(x)+Order(exp(1/x)) == Order(exp(1/x))
    assert Order(x**3)+Order(exp(2/x)) == Order(exp(2/x))
    assert Order(x**-3)+Order(exp(2/x)) == Order(exp(2/x))

def test_simple_7():
    assert 1+O(1) == O(1)
    assert 2+O(1) == O(1)
    assert x+O(1) == O(1)
    assert 1/x+O(1) == 1/x+O(1)

def test_contains_0():
    assert Order(1,x).contains(Order(1,x))
    assert Order(1,x).contains(Order(1))
    assert Order(1).contains(Order(1,x))

def test_contains_1():
    assert Order(x).contains(Order(x))
    assert Order(x).contains(Order(x**2))
    assert not Order(x**2).contains(Order(x))
    assert not Order(x).contains(Order(1/x))
    assert not Order(1/x).contains(Order(exp(1/x)))
    assert not Order(x).contains(Order(exp(1/x)))
    assert Order(1/x).contains(Order(x))
    assert Order(exp(1/x)).contains(Order(x))
    assert Order(exp(1/x)).contains(Order(1/x))
    assert Order(exp(1/x)).contains(Order(exp(1/x)))
    assert Order(exp(2/x)).contains(Order(exp(1/x)))
    assert not Order(exp(1/x)).contains(Order(exp(2/x)))

def test_contains_2():
    assert Order(x).contains(Order(y)) is None
    assert Order(x).contains(Order(y*x))
    assert Order(y*x).contains(Order(x))
    assert Order(y).contains(Order(x*y))
    assert Order(x).contains(Order(y**2*x))

def test_contains_3():
    assert Order(x*y**2).contains(Order(x**2*y)) is None
    assert Order(x**2*y).contains(Order(x*y**2)) is None

def test_add_1():
    assert Order(x+x) == Order(x)
    assert Order(3*x-2*x**2) == Order(x)
    assert Order(1+x) == Order(1,x)
    assert Order(1+1/x) == Order(1/x)
    assert Order(ln(x)+1/ln(x)) == Order(ln(x))
    assert Order(exp(1/x)+x) == Order(exp(1/x))
    assert Order(exp(1/x)+1/x**20) == Order(exp(1/x))

def test_ln_args():
    assert O(log(x)) + O(log(2*x)) == O(log(x))
    assert O(log(x)) + O(log(x**3)) == O(log(x))
    assert O(log(x*y)) + O(log(x)+log(y)) == O(log(x*y))

def test_multivar_0():
    assert Order(x*y).expr == x*y
    assert Order(x*y**2).expr == x*y**2
    assert Order(x*y,x).expr == x
    assert Order(x*y**2,y).expr == y**2
    assert Order(x*y*z).expr == x*y*z
    assert Order(x/y).expr == x/y
    assert Order(x*exp(1/y)).expr == x*exp(1/y)
    assert Order(exp(x)*exp(1/y)).expr == exp(1/y)

def test_multivar_0a():
    assert Order(exp(1/x)*exp(1/y)).expr == exp(1/x)*exp(1/y)

def test_multivar_1():
    assert Order(x+y).expr == x+y
    assert Order(x+2*y).expr == x+y
    assert (Order(x+y)+x).expr == (x+y)
    assert (Order(x+y)+x**2) == Order(x+y)
    assert (Order(x+y)+1/x) == 1/x+Order(x+y)
    assert Order(x**2+y*x).expr == x**2+y*x

def test_multivar_2():
    assert Order(x**2*y+y**2*x,x,y).expr == x**2*y+y**2*x

def test_multivar_mul_1():
    assert Order(x+y)*x == Order(x**2+y*x,x,y)

def test_multivar_3():
    assert (Order(x)+Order(y)).args in [
            (Order(x), Order(y)),
            (Order(y), Order(x))]
    assert Order(x)+Order(y)+Order(x+y) == Order(x+y)
    assert (Order(x**2*y)+Order(y**2*x)).args in [
            (Order(x*y**2), Order(y*x**2)),
            (Order(y*x**2), Order(x*y**2))]
    assert (Order(x**2*y)+Order(y*x)) == Order(x*y)

def test_issue369():
    x = Symbol('x')
    y = Symbol('y', negative=True)
    z = Symbol('z', complex=True)

    # check that Order does not modify assumptions about symbols
    Order(x)
    Order(y)
    Order(z)

    assert x.is_positive == None
    assert y.is_positive == False
    assert z.is_positive == None

    assert x.is_infinitesimal == None
    assert y.is_infinitesimal == None
    assert z.is_infinitesimal == None

def test_leading_order():
    assert (x+1+1/x**5).extract_leading_order(x) == ((1/x**5, O(1/x**5)),)
    assert (1+1/x).extract_leading_order(x) == ((1/x, O(1/x)),)
    assert (1+x).extract_leading_order(x) == ((1, O(1, x)),)
    assert (1+x**2).extract_leading_order(x) == ((1, O(1, x)),)
    assert (2+x**2).extract_leading_order(x) == ((2, O(1, x)),)
    assert (x+x**2).extract_leading_order(x) == ((x, O(x)),)

def test_leading_order2():
    assert set((2+pi+x**2).extract_leading_order(x)) == set(((pi, O(1, x)),
            (S(2), O(1, x))))
    assert set((2*x+pi*x+x**2).extract_leading_order(x)) == set(((2*x, O(x)),
            (x*pi, O(x))))

def test_order_leadterm():
    assert O(x**2)._eval_as_leading_term(x) == O(x**2)

def test_nan():
    assert not O(x).contains(nan)

def test_O1():
    assert O(1) == O(1, x)
    assert O(1) == O(1, y)
    assert hash(O(1)) == hash(O(1, x))
    assert hash(O(1)) == hash(O(1, y))
