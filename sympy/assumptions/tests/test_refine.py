from sympy.assumptions.ask import Q
from sympy.assumptions.refine import refine
from sympy.core.expr import Expr
from sympy.core.numbers import (I, Rational, nan, pi)
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.functions.elementary.complexes import (Abs, arg, im, re, sign)
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (atan, atan2, cos, sin)
from sympy.abc import w, x, y, z
from sympy.core.relational import Eq, Ne
from sympy.functions.elementary.piecewise import Piecewise
from sympy.matrices.expressions.matexpr import MatrixSymbol


def test_Abs():
    assert refine(Abs(x), Q.positive(x)) == x
    assert refine(1 + Abs(x), Q.positive(x)) == 1 + x
    assert refine(Abs(x), Q.negative(x)) == -x
    assert refine(1 + Abs(x), Q.negative(x)) == 1 - x

    assert refine(Abs(x**2)) != x**2
    assert refine(Abs(x**2), Q.real(x)) == x**2


def test_pow1():
    assert refine((-1)**x, Q.even(x)) == 1
    assert refine((-1)**x, Q.odd(x)) == -1
    assert refine((-2)**x, Q.even(x)) == 2**x

    # nested powers
    assert refine(sqrt(x**2)) != Abs(x)
    assert refine(sqrt(x**2), Q.complex(x)) != Abs(x)
    assert refine(sqrt(x**2), Q.real(x)) == Abs(x)
    assert refine(sqrt(x**2), Q.positive(x)) == x
    assert refine((x**3)**Rational(1, 3)) != x

    assert refine((x**3)**Rational(1, 3), Q.real(x)) != x
    assert refine((x**3)**Rational(1, 3), Q.positive(x)) == x

    assert refine(sqrt(1/x), Q.real(x)) != 1/sqrt(x)
    assert refine(sqrt(1/x), Q.positive(x)) == 1/sqrt(x)

    # powers of (-1)
    assert refine((-1)**(x + y), Q.even(x)) == (-1)**y
    assert refine((-1)**(x + y + z), Q.odd(x) & Q.odd(z)) == (-1)**y
    assert refine((-1)**(x + y + 1), Q.odd(x)) == (-1)**y
    assert refine((-1)**(x + y + 2), Q.odd(x)) == (-1)**(y + 1)
    assert refine((-1)**(x + 3)) == (-1)**(x + 1)

    # continuation
    assert refine((-1)**((-1)**x/2 - S.Half), Q.integer(x)) == (-1)**x
    assert refine((-1)**((-1)**x/2 + S.Half), Q.integer(x)) == (-1)**(x + 1)
    assert refine((-1)**((-1)**x/2 + 5*S.Half), Q.integer(x)) == (-1)**(x + 1)


def test_pow2():
    assert refine((-1)**((-1)**x/2 - 7*S.Half), Q.integer(x)) == (-1)**(x + 1)
    assert refine((-1)**((-1)**x/2 - 9*S.Half), Q.integer(x)) == (-1)**x

    # powers of Abs
    assert refine(Abs(x)**2, Q.real(x)) == x**2
    assert refine(Abs(x)**3, Q.real(x)) == Abs(x)**3
    assert refine(Abs(x)**2) == Abs(x)**2


def test_exp():
    x = Symbol('x', integer=True)
    assert refine(exp(pi*I*2*x)) == 1
    assert refine(exp(pi*I*2*(x + S.Half))) == -1
    assert refine(exp(pi*I*2*(x + Rational(1, 4)))) == I
    assert refine(exp(pi*I*2*(x + Rational(3, 4)))) == -I


def test_Piecewise():
    assert refine(Piecewise((1, x < 0), (3, True)), (x < 0)) == 1
    assert refine(Piecewise((1, x < 0), (3, True)), ~(x < 0)) == 3
    assert refine(Piecewise((1, x < 0), (3, True)), (y < 0)) == \
        Piecewise((1, x < 0), (3, True))
    assert refine(Piecewise((1, x > 0), (3, True)), (x > 0)) == 1
    assert refine(Piecewise((1, x > 0), (3, True)), ~(x > 0)) == 3
    assert refine(Piecewise((1, x > 0), (3, True)), (y > 0)) == \
        Piecewise((1, x > 0), (3, True))
    assert refine(Piecewise((1, x <= 0), (3, True)), (x <= 0)) == 1
    assert refine(Piecewise((1, x <= 0), (3, True)), ~(x <= 0)) == 3
    assert refine(Piecewise((1, x <= 0), (3, True)), (y <= 0)) == \
        Piecewise((1, x <= 0), (3, True))
    assert refine(Piecewise((1, x >= 0), (3, True)), (x >= 0)) == 1
    assert refine(Piecewise((1, x >= 0), (3, True)), ~(x >= 0)) == 3
    assert refine(Piecewise((1, x >= 0), (3, True)), (y >= 0)) == \
        Piecewise((1, x >= 0), (3, True))
    assert refine(Piecewise((1, Eq(x, 0)), (3, True)), (Eq(x, 0)))\
        == 1
    assert refine(Piecewise((1, Eq(x, 0)), (3, True)), (Eq(0, x)))\
        == 1
    assert refine(Piecewise((1, Eq(x, 0)), (3, True)), ~(Eq(x, 0)))\
        == 3
    assert refine(Piecewise((1, Eq(x, 0)), (3, True)), ~(Eq(0, x)))\
        == 3
    assert refine(Piecewise((1, Eq(x, 0)), (3, True)), (Eq(y, 0)))\
        == Piecewise((1, Eq(x, 0)), (3, True))
    assert refine(Piecewise((1, Ne(x, 0)), (3, True)), (Ne(x, 0)))\
        == 1
    assert refine(Piecewise((1, Ne(x, 0)), (3, True)), ~(Ne(x, 0)))\
        == 3
    assert refine(Piecewise((1, Ne(x, 0)), (3, True)), (Ne(y, 0)))\
        == Piecewise((1, Ne(x, 0)), (3, True))


def test_atan2():
    assert refine(atan2(y, x), Q.real(y) & Q.positive(x)) == atan(y/x)
    assert refine(atan2(y, x), Q.negative(y) & Q.positive(x)) == atan(y/x)
    assert refine(atan2(y, x), Q.negative(y) & Q.negative(x)) == atan(y/x) - pi
    assert refine(atan2(y, x), Q.positive(y) & Q.negative(x)) == atan(y/x) + pi
    assert refine(atan2(y, x), Q.zero(y) & Q.negative(x)) == pi
    assert refine(atan2(y, x), Q.positive(y) & Q.zero(x)) == pi/2
    assert refine(atan2(y, x), Q.negative(y) & Q.zero(x)) == -pi/2
    assert refine(atan2(y, x), Q.zero(y) & Q.zero(x)) is nan


def test_re():
    assert refine(re(x), Q.real(x)) == x
    assert refine(re(x), Q.imaginary(x)) is S.Zero
    assert refine(re(x+y), Q.real(x) & Q.real(y)) == x + y
    assert refine(re(x+y), Q.real(x) & Q.imaginary(y)) == x
    assert refine(re(x*y), Q.real(x) & Q.real(y)) == x * y
    assert refine(re(x*y), Q.real(x) & Q.imaginary(y)) == 0
    assert refine(re(x*y*z), Q.real(x) & Q.real(y) & Q.real(z)) == x * y * z


def test_im():
    assert refine(im(x), Q.imaginary(x)) == -I*x
    assert refine(im(x), Q.real(x)) is S.Zero
    assert refine(im(x+y), Q.imaginary(x) & Q.imaginary(y)) == -I*x - I*y
    assert refine(im(x+y), Q.real(x) & Q.imaginary(y)) == -I*y
    assert refine(im(x*y), Q.imaginary(x) & Q.real(y)) == -I*x*y
    assert refine(im(x*y), Q.imaginary(x) & Q.imaginary(y)) == 0
    assert refine(im(1/x), Q.imaginary(x)) == -I/x
    assert refine(im(x*y*z), Q.imaginary(x) & Q.imaginary(y)
        & Q.imaginary(z)) == -I*x*y*z


def test_complex():
    assert refine(re(1/(x + I*y)), Q.real(x) & Q.real(y)) == \
        x/(x**2 + y**2)
    assert refine(im(1/(x + I*y)), Q.real(x) & Q.real(y)) == \
        -y/(x**2 + y**2)
    assert refine(re((w + I*x) * (y + I*z)), Q.real(w) & Q.real(x) & Q.real(y)
        & Q.real(z)) == w*y - x*z
    assert refine(im((w + I*x) * (y + I*z)), Q.real(w) & Q.real(x) & Q.real(y)
        & Q.real(z)) == w*z + x*y


def test_sign():
    x = Symbol('x', real = True)
    assert refine(sign(x), Q.positive(x)) == 1
    assert refine(sign(x), Q.negative(x)) == -1
    assert refine(sign(x), Q.zero(x)) == 0
    assert refine(sign(x), True) == sign(x)
    assert refine(sign(Abs(x)), Q.nonzero(x)) == 1

    x = Symbol('x', imaginary=True)
    assert refine(sign(x), Q.positive(im(x))) == S.ImaginaryUnit
    assert refine(sign(x), Q.negative(im(x))) == -S.ImaginaryUnit
    assert refine(sign(x), True) == sign(x)

    x = Symbol('x', complex=True)
    assert refine(sign(x), Q.zero(x)) == 0

def test_arg():
    x = Symbol('x', complex = True)
    assert refine(arg(x), Q.positive(x)) == 0
    assert refine(arg(x), Q.negative(x)) == pi

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

def test_issue_refine_9384():
    assert refine(Piecewise((1, x < 0), (0, True)), Q.positive(x)) == 0
    assert refine(Piecewise((1, x < 0), (0, True)), Q.negative(x)) == 1
    assert refine(Piecewise((1, x > 0), (0, True)), Q.positive(x)) == 1
    assert refine(Piecewise((1, x > 0), (0, True)), Q.negative(x)) == 0


def test_eval_refine():
    class MockExpr(Expr):
        def _eval_refine(self, assumptions):
            return True

    mock_obj = MockExpr()
    assert refine(mock_obj)

def test_refine_issue_12724():
    expr1 = refine(Abs(x * y), Q.positive(x))
    expr2 = refine(Abs(x * y * z), Q.positive(x))
    assert expr1 == x * Abs(y)
    assert expr2 == x * Abs(y * z)
    y1 = Symbol('y1', real = True)
    expr3 = refine(Abs(x * y1**2 * z), Q.positive(x))
    assert expr3 == x * y1**2 * Abs(z)


def test_matrixelement():
    x = MatrixSymbol('x', 3, 3)
    i = Symbol('i', positive = True)
    j = Symbol('j', positive = True)
    assert refine(x[0, 1], Q.symmetric(x)) == x[0, 1]
    assert refine(x[1, 0], Q.symmetric(x)) == x[0, 1]
    assert refine(x[i, j], Q.symmetric(x)) == x[j, i]
    assert refine(x[j, i], Q.symmetric(x)) == x[j, i]


def test_sin_cos():
    n = Symbol('n')
    assert refine(cos(n*pi/2), Q.odd(n)) == 0
    assert refine(cos(n*pi), Q.even(n)) == 1
    assert refine(cos(n*pi), Q.odd(n)) == -1
    assert refine(sin(n*pi), Q.integer(n)) == 0
    assert refine(sin(n*pi/2), Q.odd(n) & Q.even((n-1)/2)) == 1
    assert refine(sin(n*pi/2), Q.odd(n) & Q.odd((n-1)/2)) == -1
    assert refine(cos(n*pi), Q.integer(n)) == (-1)**n
    assert refine(sin(n*pi/2), Q.even(n)) == 0
    assert refine(cos(n*pi/2), Q.even(n)) == (-1)**(n/2)
    assert refine(sin(n*pi/2), Q.odd(n)) == (-1)**((n + 3)/2)
    assert refine(cos(n*pi/2), Q.odd(n)) == 0
    assert refine(sin(x + n*pi), Q.integer(n)) == ((-1)**n) * sin(x)
    assert refine(cos(x + n*pi), Q.integer(n)) == ((-1)**n) * cos(x)
    assert refine(sin(x + n*pi), Q.even(n)) == sin(x)
    assert refine(cos(x + n*pi), Q.even(n)) == cos(x)
    assert refine(sin(x + n*pi), Q.odd(n)) == -sin(x)
    assert refine(cos(x + n*pi), Q.odd(n)) == -cos(x)
    assert refine(sin(x - n*pi), Q.odd(n)) == -sin(x)
    assert refine(cos(x - n*pi), Q.even(n)) == cos(x)
    assert refine(sin(x + n*pi/2), Q.even(n)) == ((-1)**(n/2)) * sin(x)
    assert refine(cos(x + n*pi/2), Q.even(n)) == ((-1)**(n/2)) * cos(x)
    assert refine(sin(x + n*pi/2), Q.odd(n)) == ((-1)**((n + 3)/2)) * cos(x)
    assert refine(cos(x + n*pi/2), Q.odd(n)) == ((-1)**((n + 1)/2)) * sin(x)
    assert refine(sin(x - n*pi/2), Q.odd(n)) == ((-1)**((n + 3)/2)) * -cos(x)
    assert refine(cos(x - n*pi / 2), Q.even(n)) == ((-1)**(n/2)) * cos(x)
    assert refine(sin(x + y + 2*n*pi), Q.integer(n)) == sin(x + y)
    assert refine(cos(x + y + 2*n*pi), Q.integer(n)) == cos(x + y)
    assert refine(sin(x + n*pi), Q.zero(n)) == sin(x)
    assert refine(sin(x + n*pi), Q.zero(-n)) == sin(x)
    assert refine(cos(x + n*pi/2), Q.integer(n)) == cos(x + n*pi/2)
    assert refine(cos(x + y + n*pi/2), Q.integer(n)) == cos(x + y + n*pi/2)
    m = Symbol('m')
    assert refine(cos(x + n*pi + m*pi / 2), Q.integer(n) & Q.even(m)) == \
        (-1)**(n + m / 2) * cos(x)
    assert refine(cos(x + n*pi + m*pi / 2), Q.integer(n) & Q.odd(m)) == \
        (-1)**(n + (m + 1)/2) * sin(x)
    assert refine(cos(x + n*pi + m*pi / 2), Q.integer(n) & Q.integer(m)) == \
        (-1)**(n) * cos(x + m*pi / 2)
    assert refine(cos(x + (2*n + 1)*pi + m*pi / 2), \
        Q.integer(n) & Q.integer(m)) == \
        - cos(x + m*pi / 2)
    assert refine(sin(x - (2*n)*pi + m*pi/2), \
        Q.integer(n) & Q.integer(m)) == \
        sin(x + m*pi / 2)
    k = Symbol('k')
    assert refine(cos(x + n*pi + k*pi/2 + m*pi/2), \
                  Q.integer(n) & Q.odd(k) & Q.integer(m)) == \
        (-1)**(n + (k + 1)/2) * sin(x + m*pi/2)
    assert refine(sin(x + n*pi + k*pi/2 + m*pi/2), \
                  Q.integer(n) & Q.odd(k) & Q.integer(m)) == \
        (-1)**(n + (k + 3)/2) * cos(x + m*pi/2)
    assert refine(cos(x + n*pi/2 + k*pi/2 + m*pi/2), \
                  Q.odd(n) & Q.odd(k) & Q.integer(m)) == \
        (-1)**((n + k)/2) * cos(x + m*pi/2)


def test_log():
    # Power rule: log(x**p) -> p*log(x) when x > 0 and p is real
    assert refine(log(x**2), Q.positive(x)) == 2*log(x)
    assert refine(log(x**y), Q.positive(x) & Q.real(y)) == y*log(x)
    assert refine(log(x**3), Q.positive(x)) == 3*log(x)

    # Power rule must NOT apply with imaginary exponents
    assert refine(log(x**(2*I)), Q.positive(x)) == log(x**(2*I))
    assert refine(log(x**y), Q.positive(x)) == log(x**y)

    # Mixed exponent: split real parts from non-real parts
    assert (refine(log(x**(2*I + y)), Q.positive(x) & Q.real(y)) ==
        y*log(x) + log(x**(2*I)))
    # Also works with Q.positive(y) since positive implies real
    assert (refine(log(x**(2*I + y)), Q.positive(x) & Q.positive(y)) ==
        y*log(x) + log(x**(2*I)))

    # Product rule: log(x*y) -> log(x) + log(y) when both positive
    assert refine(log(x*y), Q.positive(x) & Q.positive(y)) == log(x) + log(y)
    assert (refine(log(x*y*z), Q.positive(x) & Q.positive(y) & Q.positive(z)) ==
        log(x) + log(y) + log(z))

    # Quotient rule: log(x/y) -> log(x) - log(y) when both positive
    # (x/y is represented as x * y**(-1))
    assert refine(log(x/y), Q.positive(x) & Q.positive(y)) == log(x) - log(y)

    # Inverse: log(exp(x)) -> x when x is real
    assert refine(log(exp(x)), Q.real(x)) == x
    assert refine(log(exp(y)), Q.real(y)) == y

    # No simplification without proper assumptions
    assert refine(log(x**2), True) == log(x**2)
    assert refine(log(x*y), Q.positive(x)) == log(x) + log(y)  # partial split
    assert refine(log(exp(x)), True) == log(exp(x))  # need Q.real

    # Combined product + power: log(x**2 * y**3) with both positive
    assert (refine(log(x**2 * y**3), Q.positive(x) & Q.positive(y)) ==
        2*log(x) + 3*log(y))

    # Nested: log(exp(x)) inside a product
    assert (refine(log(x * exp(y)), Q.positive(x) & Q.real(y)) ==
        log(x) + y)

    # log(1/x) with positive x -> -log(x)
    assert refine(log(S.One/x), Q.positive(x)) == -log(x)

    # Negative base should NOT simplify power rule
    assert refine(log(x**2), Q.negative(x)) == log(x**2)

    # Both terms negative should NOT split
    assert refine(log(x*y), Q.negative(x) & Q.negative(y)) == log(x*y)
    
    # One negative, one positive — positive term still extracts
    assert refine(log(x*y), Q.positive(x) & Q.negative(y)) == log(x) + log(y)
    
    # Complex exponent with no real part should NOT split
    assert refine(log(x**(3*I)), Q.positive(x)) == log(x**(3*I))
    assert refine(log(x**(I*y)), Q.positive(x) & Q.real(y)) == log(x**(I*y))
    
    # exp with complex argument should NOT fully simplify
    assert refine(log(exp(I*x)), Q.real(x)) == log(exp(I*x))
    assert refine(log(exp(x + I)), Q.real(x)) == x + log(exp(I))
    
    # Unknown sign should NOT split product
    assert refine(log(x*y), True) == log(x*y)
    assert refine(log(x*y), Q.real(x) & Q.real(y)) == log(x*y)
    
    # Zero should NOT simplify (undefined)
    assert refine(log(0), True) == log(0)
    
    # (-x)**2 is canonicalized to x**2 by SymPy, so it simplifies
    assert refine(log((-x)**2), Q.positive(x)) == 2*log(x)
    
    # Mixed known/unknown in product
    assert refine(log(x*y), Q.positive(x)) == log(x) + log(y)  # partial OK
    assert refine(log(x*y), Q.negative(x)) == log(x*y)  # no split
    
    # Non-real exponent without Q.real assumption
    assert refine(log(x**y), Q.positive(x)) == log(x**y)  # y could be complex
    
    # Multiple imaginary terms in exponent
    assert refine(log(x**(I + 2*I)), Q.positive(x)) == log(x**(I + 2*I))
    
    # log(exp(x)) without Q.real should NOT simplify
    assert refine(log(exp(x)), True) == log(exp(x))
    assert refine(log(exp(x)), Q.positive(x)) == x  # positive implies real
    
    # Nested log should NOT simplify
    assert refine(log(log(x)), Q.positive(x)) == log(log(x))
    
    # Sum inside log (not product) should NOT split
    assert refine(log(x + y), Q.positive(x) & Q.positive(y)) == log(x + y)
    
    # Power with negative exponent
    assert refine(log(x**(-2)), Q.positive(x)) == -2*log(x)
    assert refine(log(x**(-y)), Q.positive(x) & Q.real(y)) == -y*log(x)
    
    # Four terms with mixed assumptions
    assert (refine(log(x*y*z*w), Q.positive(x) & Q.negative(y)) ==
        log(x) + log(y*z*w))  # Only x extracts
    
    # All four positive should fully split
    assert (refine(log(x*y*z*w), Q.positive(x) & Q.positive(y) & 
                   Q.positive(z) & Q.positive(w)) ==
        log(x) + log(y) + log(z) + log(w))
    
    # Three positive, one unknown
    assert (refine(log(x*y*z*w), Q.positive(x) & Q.positive(y) & 
                   Q.positive(z)) ==
        log(x) + log(y) + log(z) + log(w))  # w stays inside
    
    # exp with multiple real variables
    assert refine(log(exp(x + y)), Q.real(x) & Q.real(y)) == x + y
    assert refine(log(exp(x + y + z)), Q.real(x) & Q.real(y) & Q.real(z)) == x + y + z
    
    # exp with mixed real/imaginary
    assert (refine(log(exp(x + y + I*z)), Q.real(x) & Q.real(y) & Q.real(z)) ==
        x + y + log(exp(I*z)))
    
    # Power with multiple variables in exponent
    assert (refine(log(x**(y + z)), Q.positive(x) & Q.real(y) & Q.real(z)) ==
        (y + z)*log(x))
    assert (refine(log(x**(y + I*z)), Q.positive(x) & Q.real(y) & Q.real(z)) ==
        y*log(x) + log(x**(I*z)))

    # Power rule needs Q.positive(x) AND Q.real(exponent)
    # -- omit Q.positive(x): only Q.real
    assert refine(log(x**2), Q.real(x)) == log(x**2)
    assert refine(log(x**y), Q.real(y)) == log(x**y)
    assert refine(log(x**y), Q.real(x) & Q.real(y)) == log(x**y)

    # Inverse rule needs Q.real(x)
    # -- wrong assumption: imaginary or complex instead of real
    assert refine(log(exp(x)), Q.imaginary(x)) == log(exp(x))
    assert refine(log(exp(x)), Q.complex(x)) == log(exp(x))

    # Mixed exponent needs Q.positive(x) AND Q.real(y)
    # -- omit Q.real(y)
    assert refine(log(x**(y + 2*I)), Q.positive(x)) == log(x**(y + 2*I))
    # -- omit Q.positive(x)
    assert refine(log(x**(y + 2*I)), Q.real(y)) == log(x**(y + 2*I))
    # -- omit both
    assert refine(log(x**(y + 2*I)), True) == log(x**(y + 2*I))

    # exp with mixed args needs Q.real for the real parts
    # -- no assumptions
    assert refine(log(exp(x + I)), True) == log(exp(x + I))
    # -- wrong variable assumed real
    assert refine(log(exp(x + I*y)), Q.real(y)) == log(exp(x + I*y))

    # Quotient rule needs both Q.positive
    # -- only Q.positive(y), omit Q.positive(x)
    assert refine(log(x/y), Q.positive(y)) == log(x) - log(y)
    # -- no assumptions
    assert refine(log(x/y), True) == log(x/y)

    # Multi-term product: no assumptions
    assert refine(log(x*y*z), True) == log(x*y*z)

    # exp multi real: omit one Q.real
    assert refine(log(exp(x + y)), Q.real(x)) == x + log(exp(y))
    # -- no assumptions
    assert refine(log(exp(x + y)), True) == log(exp(x + y))

    # exp mixed real/imag: omit one Q.real
    assert (refine(log(exp(x + y + I*z)), Q.real(x)) ==
        x + log(exp(y + I*z)))
    # -- no assumptions
    assert refine(log(exp(x + y + I*z)), True) == log(exp(x + y + I*z))

    # Power multi var exponent: omit assumptions
    # -- omit Q.real(y) and Q.real(z)
    assert refine(log(x**(y + z)), Q.positive(x)) == log(x**(y + z))
    # -- omit Q.positive(x)
    assert refine(log(x**(y + z)), Q.real(y) & Q.real(z)) == log(x**(y + z))
    # -- omit Q.real(z) from mixed
    assert refine(log(x**(y + I*z)), Q.positive(x)) == log(x**(y + I*z))