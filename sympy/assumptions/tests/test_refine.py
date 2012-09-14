from sympy import Abs, exp, Expr, I, pi, Q, Rational, refine, S, sqrt
from sympy.abc import x, y, z

def test_Abs():
    assert refine(Abs(x), Q.positive(x)) == x
    assert refine(1+Abs(x), Q.positive(x)) == 1+x
    assert refine(Abs(x), Q.negative(x)) == -x
    assert refine(1+Abs(x), Q.negative(x)) == 1-x

    assert refine(Abs(x**2)) != x**2
    assert refine(Abs(x**2), Q.real(x)) == x**2

def test_pow():
    assert refine((-1)**x, Q.even(x)) == 1
    assert refine((-1)**x, Q.odd(x)) == -1
    assert refine((-2)**x, Q.even(x)) == 2**x

    # nested powers
    assert refine(sqrt(x**2)) != Abs(x)
    assert refine(sqrt(x**2), Q.complex(x)) != Abs(x)
    assert refine(sqrt(x**2), Q.real(x)) == Abs(x)
    assert refine(sqrt(x**2), Q.positive(x)) == x
    assert refine((x**3)**(S(1)/3)) != x

    assert refine((x**3)**(S(1)/3), Q.real(x)) != x
    assert refine((x**3)**(S(1)/3), Q.positive(x)) == x

    assert refine(sqrt(1/x), Q.real(x)) != 1/sqrt(x)
    assert refine(sqrt(1/x), Q.positive(x)) == 1/sqrt(x)

    # powers of (-1)
    assert refine((-1)**(x+y), Q.even(x)) == (-1)**y
    assert refine((-1)**(x+y+z), Q.odd(x) & Q.odd(z))==(-1)**y
    assert refine((-1)**(x+y+1), Q.odd(x))==(-1)**y
    assert refine((-1)**(x+y+2), Q.odd(x))==(-1)**(y+1)
    assert refine((-1)**(x+3)) == (-1)**(x+1)


def test_exp():
    assert refine(exp(pi*I*2*x), Q.integer(x)) == 1
    assert refine(exp(pi*I*2*(x+Rational(1,2))), Q.integer(x)) == -1
    assert refine(exp(pi*I*2*(x+Rational(1,4))), Q.integer(x)) == I
    assert refine(exp(pi*I*2*(x+Rational(3,4))), Q.integer(x)) == -I


def test_func_args():
    class MyClass(Expr):
        # A class with nontrivial .func

        def __init__(self, *args):
            self.my_member = ""

        @property
        def func(self):
            def my_func(*args):
                obj = MyClass(*args)
                obj.my_member = self.my_member
                return obj
            return my_func

    x = MyClass()
    x.my_member = "A very important value"
    assert x.my_member == refine(x).my_member
