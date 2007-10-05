from sympy import *
from sympy.printing.pretty import pretty
from sympy.utilities.pytest import XFAIL

x = Symbol('x')
y = Symbol('y')

def test_pretty_basic():
    # Simple numbers/symbols
    assert pretty( -Rational(1)/2 ) == '-1/2'
    assert pretty( -Rational(13)/22 ) == '  13\n- --\n  22'
    assert pretty( oo ) == 'oo'

    # Powers
    assert pretty( (x**2) ) == ' 2\nx '
    assert pretty( 1/x ) == '1\n-\nx'

    # Sums of terms
    assert pretty( (x**2 + x + 1))  == '         2\n1 + x + x '
    assert pretty( 1-x ) == '1 - x'
    assert pretty( 1-2*x ) == '1 - 2*x'
    assert pretty( 1-Rational(3,2)*y/x ) == '    3*y\n1 - ---\n    2*x'

    # Multiplication
    assert pretty( x/y ) == 'x\n-\ny'
    assert pretty( -x/y ) == '-x\n--\ny '
    assert pretty( (x+2)/y ) == '2 + x\n-----\n  y  '
    assert pretty( (1+x)*y ) in ['(1 + x)*y', 'y*(1 + x)']

    # Check for proper placement of negative sign
    assert pretty( -5*x/(x+10) ) == ' -5*x \n------\n10 + x'
    assert pretty( 1 - Rational(3,2)*(x+1) ) == '       3*x\n-1/2 - ---\n        2 '

def test_pretty_relational():
    assert pretty(x == y) == 'x = y'
    assert pretty(x <= y) == 'x <= y'
    assert pretty(x > y) == 'y < x'
    assert pretty(x/(y+1) != y**2) == '  x       2\n----- != y \n1 + y      '

def test_pretty_unicode():
    assert pretty( oo, True ) == u'\u221e'
    assert pretty( pi, True ) == u'\u03c0'
    assert pretty( pi+2*x, True ) == u'\u03c0 + 2*x'
    assert pretty( pi**2+exp(x), True ) == u' 2    x\n\u03c0  + \u212f '
    assert pretty( x != y, True ) == u'x \u2260 y'

def test_pretty_unicode_defaults():
    use_unicode = pprint_use_unicode(True)
    assert pretty(Symbol('alpha')) == u'\u03b1'
    pprint_use_unicode(False)
    assert pretty(Symbol('alpha')) == 'alpha'

    pprint_use_unicode(use_unicode)


def test_pretty_functions():
    # Simple
    assert pretty( (2*x + exp(x)) ) in [' x      \ne  + 2*x', '       x\n2*x + e ']
    assert pretty( sqrt(2) ) == '  ___\n\\/ 2 '
    assert pretty( sqrt(2+pi) ) == '  ________\n\\/ 2 + pi '
    assert pretty(abs(x)) == '|x|'
    assert pretty(abs(x/(x**2+1))) == '|  x   |\n|------|\n|     2|\n|1 + x |'

    # Univariate/Multivariate functions
    f = Function('f')
    assert pretty(f(x)) == 'f(x)'
    assert pretty(f(x, y)) == 'f(x, y)'
    assert pretty(f(x/(y+1), y)) == '    x      \nf(-----, y)\n  1 + y    '

    # Nesting of square roots
    assert pretty( sqrt((sqrt(x+1))+1) ) == '   _______________\n  /       _______ \n\\/  1 + \\/ 1 + x  '
    # Function powers
    assert pretty( sin(x)**2 ) == '   2   \nsin (x)'

    # Conjugates
    a,b = map(Symbol, 'ab')
    #assert pretty( conjugate(a+b*I) ) == '_     _\na - I*b'
    #assert pretty( conjugate(exp(a+b*I)) ) == ' _     _\n a - I*b\ne       '

def test_pretty_derivatives():
    # Simple
    f_1 = Derivative(log(x), x, evaluate=False)
    assert pretty(f_1) == 'd         \n--(log(x))\ndx        '

    f_2 = Derivative(log(x), x, evaluate=False) + x
    assert pretty(f_2) == '    d         \nx + --(log(x))\n    dx        '

    # Multiple symbols
    f_3 = Derivative(log(x) + x**2, x, y, evaluate=False)
    assert pretty(f_3) == '   2              \n  d  / 2         \\\n-----\\x  + log(x)/\ndy dx             '

    f_4 = Derivative(2*x*y, y, x, evaluate=False) + x**2
    assert pretty(f_4) == '        2        \n 2     d         \nx  + -----(2*x*y)\n     dx dy       '

def test_pretty_integrals():
    # Simple
    f_1 = Integral(log(x), x)
    assert pretty(f_1) == '  /         \n |          \n | log(x) dx\n |          \n/           '

    f_2 = Integral(x**2, x)
    assert pretty(f_2) == '  /     \n |      \n |  2   \n | x  dx\n |      \n/       '
    # Double nesting of pow
    f_3 = Integral(x**(2**x), x)
    assert pretty(f_3) == '  /        \n |         \n |  / x\\   \n |  \\2 /   \n | x     dx\n |         \n/          '

    # Definite integrals
    f_4 = Integral(x**2, (x,1,2))
    assert pretty(f_4) == '  2      \n  /      \n |       \n |   2   \n |  x  dx\n |       \n/        \n1        '

    f_5 = Integral(x**2, (x,Rational(1,2),10))
    assert pretty(f_5) == ' 10      \n  /      \n |       \n |   2   \n |  x  dx\n |       \n/        \n1/2      '

    # Nested integrals
    f_6 = Integral(x**2*y**2, x,y)
    assert pretty(f_6) == '  /  /           \n |  |            \n |  |  2  2      \n |  | x *y  dx dy\n |  |            \n/  /             '

@XFAIL
def test_pretty_limits():
    assert pretty( limit(x, x, oo, evaluate=False) ) == ' lim x\nx->oo '
    assert pretty( limit(x**2, x, 0, evaluate=False) ) == '     2\nlim x \nx->0  '  
