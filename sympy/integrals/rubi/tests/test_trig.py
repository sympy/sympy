from random import randint
from sympy import S, Symbol, sin, cos
from sympy.integrals.rubi.SineIntegrationRules import intsin5, intsin12

# The various functions in SineIntegrationRules return the antiderivatives of
# the following expressions:
#   intsin5(a,b,c,n,x)
#       (c*sin(a+b*x))**n
#   intsin6(a,b,c,d,n,x)
#       (a+b*sin(c+d*x))**n
#   intsin9(a,b,c,d,e,f,m,n,x)
#       (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n
#   intsin(a,b,c,d,e,f,A,B,m,n,x)
#       x*(a+b*sin(e+f*x))**m*(A+B*sin(e+f*x))*(c+d*sin(e+f*x))**n
#   intsin12(a,b,c,d,e,f,A,B,C,m,n,x)
#       (a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n*(A+B*sin(e+f*x)+C*sin(e+f*x)**2)
# where the parameters (a,b,c,d,e,f,A,B,C) and exponents (m,n) are arbitrary
# numeric or symbolic expressions, including zeros and ones.

def test_intsin5():
    x = Symbol("x")
    assert intsin5(S(1), S(2), S(3), S(2), x) == 9*x/2 - 9*sin(2*x + 1)*cos(2*x + 1)/4
    assert intsin5(S(1), S(2), S(3), -S(2), x) == -cos(2*x + 1)/(18*sin(2*x + 1))

def test_intsin12():
    x = Symbol("x")
    assert intsin12(S(0), S(3), S(1), S(0), S(1), S(2), S(1), S(0), S(0), S(2), S(0), x) == 9*x/2 - 9*sin(2*x + 1)*cos(2*x + 1)/4

def is_zero(e, n=10, prec=100):
    """
    Tests that the expression 'e' is zero.

    It is just a probabilistic test by evaluating the expression 'e' at various
    points and making sure the values are all sufficiently close to zero. If
    the function returns False, then the expressions are *not* equal. If the
    function returns True, the expressions *might* be equal.

    n ...... the number of substitutions to try
    prec ... the value of expression must be less than 10^(-prec) in order to
             be zero

    The expression can contain any number of symbols (the test is slower the
    more symbols it has). is_zero() can handle any number of symbols.
    """
    symbols = list(e.atoms(Symbol))
    for i in range(n):
        d = {}
        for j in range(len(symbols)):
            x0 = S(randint(-10000, 10000))/500
            d[symbols[j]] = x0
        if (abs(e.subs(d).n(prec)) > S(10)**(-prec)):
            return False
    return True

def test_is_zero():
    assert is_zero(S(0))
    assert not is_zero(S(1))
    assert not is_zero(S(-1))
    assert is_zero(sin(1)**2+cos(1)**2-1)
    assert not is_zero(sin(1)**2+cos(1)**2)

    x = Symbol("x")
    assert is_zero((x+1)**2-x**2-2*x-1)
    assert not is_zero((x+1)**2-x**2-2*x)
    assert is_zero(sin(x)**2+cos(x)**2-1)
    assert not is_zero(sin(x)**2+cos(x)**2)

    y = Symbol("y")
    assert is_zero(sin(x)**2+cos(x)**2-1+sin(y)**2+cos(y)**2-1)
    assert not is_zero(sin(x)**2+cos(x)**2-1+sin(y)**2+cos(y)**2)

    z = Symbol("z")
    assert is_zero(sin(x)**2+cos(x)**2-1+sin(y)**2+cos(y)**2-1+sin(z)**2+cos(z)**2-1)
    assert not is_zero(sin(x)**2+cos(x)**2-1+sin(y)**2+cos(y)**2-1+sin(z)**2+cos(z)**2)

def check_intsin12(a,b,c,d,e,f,A,B,C,m,n,x):
    expression = (a+b*sin(e+f*x))**m * (c+d*sin(e+f*x))**n * \
            (A+B*sin(e+f*x)+C*sin(e+f*x)**2)
    integral = intsin12(a,b,c,d,e,f,A,B,C,m,n,x)
    zero = expression - integral.diff(x)
    return is_zero(zero)

def test_auto():
    x = Symbol("x")
    args = [S(0), S(3), S(1), S(0), S(1), S(2), S(1), S(0), S(0), S(2), S(0), x]
    assert intsin12(*args) == 9*x/2 - 9*sin(2*x + 1)*cos(2*x + 1)/4
    assert check_intsin12(*args)
    args = [S(0), S(3), S(1), S(0), S(1), S(2), S(1), S(1), S(0), S(2), S(0), x]
    assert check_intsin12(*args)
    args = [S(0), S(3), S(1), S(0), S(1), S(2), S(1), S(1), S(1), S(2), S(0), x]
    assert check_intsin12(*args)
    args = [S(1), S(3), S(1), S(1), S(1), S(2), S(1), S(1), S(1), S(2), S(2), x]
    assert check_intsin12(*args)

    a = Symbol("a")
    args = [S(a), S(3), S(1), S(1), S(1), S(2), S(1), S(1), S(1), S(2), S(2), x]
    assert check_intsin12(*args)

    b = Symbol("b")
    args = [S(a), S(3), S(1), S(1), S(1), S(2), S(1), S(b), S(1), S(2), S(2), x]
    assert check_intsin12(*args)

    c = Symbol("c")
    args = [S(a), S(3), S(1), S(1), S(1), S(c), S(1), S(b), S(1), S(2), S(2), x]
    assert check_intsin12(*args)
