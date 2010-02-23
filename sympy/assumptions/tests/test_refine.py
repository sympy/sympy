from sympy import S, symbols, Assume, exp, pi, sqrt, Rational, I, Q, refine
from sympy.utilities.pytest import XFAIL

def test_abs():
    x = symbols('x')
    assert refine(abs(x), Assume(x, Q.positive)) == x
    assert refine(1+abs(x), Assume(x, Q.positive)) == 1+x
    assert refine(abs(x), Assume(x, Q.negative)) == -x
    assert refine(1+abs(x), Assume(x, Q.negative)) == 1-x

    assert refine(abs(x**2)) != x**2
    assert refine(abs(x**2), Assume(x, Q.real)) == x**2

def test_pow():
    x, y, z = symbols('x y z')
    assert refine((-1)**x, Assume(x, Q.even)) == 1
    assert refine((-1)**x, Assume(x, Q.odd)) == -1
    assert refine((-2)**x, Assume(x, Q.even)) == 2**x

    # nested powers
    assert refine(sqrt(x**2)) != abs(x)
    assert refine(sqrt(x**2), Assume(x, Q.complex)) != abs(x)
    assert refine(sqrt(x**2), Assume(x, Q.real)) == abs(x)
    assert refine(sqrt(x**2), Assume(x, Q.positive)) == x
    assert refine((x**3)**(S(1)/3)) != x

    assert refine((x**3)**(S(1)/3), Assume(x, Q.real)) != x
    assert refine((x**3)**(S(1)/3), Assume(x, Q.positive)) == x

    assert refine(sqrt(1/x), Assume(x, Q.real)) != 1/sqrt(x)
    assert refine(sqrt(1/x), Assume(x, Q.positive)) == 1/sqrt(x)

    # powers of (-1)
    assert refine((-1)**(x+y), Assume(x, Q.even)) == (-1)**y
    assert refine((-1)**(x+y+z), Assume(x, Q.odd)&Assume(z, Q.odd))==(-1)**y
    assert refine((-1)**(x+y+1), Assume(x, Q.odd))==(-1)**y
    assert refine((-1)**(x+y+2), Assume(x, Q.odd))==(-1)**(y+1)
    assert refine((-1)**(x+3)) == (-1)**(x+1)


def test_exp():
    x = symbols('x')
    assert refine(exp(pi*I*2*x), Assume(x, Q.integer)) == 1
    assert refine(exp(pi*I*2*(x+Rational(1,2))), Assume(x, Q.integer)) == -1
    assert refine(exp(pi*I*2*(x+Rational(1,4))), Assume(x, Q.integer)) == I
    assert refine(exp(pi*I*2*(x+Rational(3,4))), Assume(x, Q.integer)) == -I
