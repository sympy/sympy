from sympy import *
from decimal import Decimal

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')


def feq(a,b):
    t = 1E-15
    return -t < a-b < t

def test_cos():
    assert cos(0) == 1
    assert cos(pi) == -1
    assert cos(-pi) == -1
    assert cos(2*pi) == 1
    assert cos(3*pi) == -1
    assert cos(10**157*pi) == 1
    assert cos((10**157+1)*pi) == -1

    assert cos(pi/2) == 0
    assert cos((10**157+123)*pi/2) == 0

    #assert cos(pi/3) == Rational(1, 2)
    #assert cos(2*pi/3) == -Rational(1, 2)
    #assert cos(-2*pi/3) == -Rational(1, 2)
    #assert cos(4*pi/3) == -Rational(1, 2)
    #assert cos(5*pi/3) == Rational(1, 2)
    #assert cos(7*pi/3) == Rational(1, 2)
    #assert cos(8*pi/3) == -Rational(1, 2)

    #assert cos(pi/4) == Rational(1, 2)*sqrt(2)
    #assert cos(-pi/4) == Rational(1, 2)*sqrt(2)
    #assert cos(3*pi/4) == -Rational(1, 2)*sqrt(2)
    #assert cos(5*pi/4) == -Rational(1, 2)*sqrt(2)
    #assert cos(7*pi/4) == Rational(1, 2)*sqrt(2)

    #assert cos(pi/6) == Rational(1, 2)*sqrt(3)
    #assert cos(5*pi/6) == -Rational(1, 2)*sqrt(3)

    assert feq( float(cos(1)), 0.54030230586813977 )
    assert feq( float(cos(1) + cos(2)), 0.12415546932099736 )
    assert feq( float(cos(1)*cos(2)*cos(3)), 0.22259495730990297 )
    assert cos(1).evalf(precision=28) == Real("0.5403023058681397174009366074")

    assert cos(-x).expand() == cos(x)
    assert cos(-x-y).expand() == cos(x)*cos(y) - sin(x)*sin(y)
    assert cos(x+y-z).expand() == cos(x)*cos(y)*cos(z) + cos(x)*sin(y)*sin(z) - sin(x)*sin(y)*cos(z) + sin(x)*cos(y)*sin(z)

    assert cos(2*x).expand() == 2*cos(x)**2 - 1
    assert cos(-3*x).expand() == 4*cos(x)**3 - 3*cos(x)

def test_sin():
    assert sin(0) == 0
    assert sin(pi) == 0
    assert sin(-pi) == 0
    assert sin(2*pi) == 0
    assert sin(-3*10**73*pi) == 0
    assert sin(7*10**103*pi) == 0

    assert sin(pi/2) == 1
    assert sin(-pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(7*pi/2) == -1

    #assert sin(pi/3) == Rational(1, 2)*sqrt(3)
    #assert sin(-2*pi/3) == -Rational(1, 2)*sqrt(3)

    #assert sin(pi/4) == Rational(1, 2)*sqrt(2)
    #assert sin(-pi/4) == -Rational(1, 2)*sqrt(2)
    #assert sin(17*pi/4) == Rational(1, 2)*sqrt(2)
    #assert sin(-3*pi/4) == -Rational(1, 2)*sqrt(2)

    #assert sin(pi/6) == Rational(1, 2)
    #assert sin(-pi/6) == -Rational(1, 2)
    #assert sin(7*pi/6) == -Rational(1, 2)
    #assert sin(-5*pi/6) == -Rational(1, 2)

    assert feq( float(sin(1)), 0.8414709848078965 )
    assert sin(1).evalf(precision=28) == Real("0.8414709848078965066525023216")

    assert sin(-x).expand() == -sin(x)
    assert sin(-x-y).expand() == -sin(x)*cos(y) - cos(x)*sin(y)
    assert sin(x+y+z).expand() == sin(x)*cos(y)*cos(z) - sin(x)*sin(y)*sin(z) + cos(x)*sin(y)*cos(z) + cos(x)*cos(y)*sin(z)

    assert sin(2*x).expand() == 2*sin(x)*cos(x)
    assert sin(-3*x).expand() == -4*sin(x)*cos(x)**2 + sin(x)

def _test_tan():
    # XXX These all fail since expansion doesn't wotk
    assert tan(-x).expand() == -tan(x)
    assert tan(-x-y).expand() == (-tan(x) - tan(y)) / (1 - tan(x)*tan(y))
    assert tan(x+y+z).expand() == \
           (tan(x) + tan(y) + tan(z) - tan(x)*tan(y)*tan(z)) \
           / (1 - tan(x)*tan(y) - tan(x)*tan(z) - tan(y)*tan(z)) 

    assert tan(2*x).expand() == 2*tan(x) / (1 - tan(x)**2)
    assert tan(-3*x).expand() == -((3*tan(x) - tan(x)**3) / (1 - 3*tan(x)**2))

def _test_tan_const():
    assert tan(0)       == 0
    # XXX This fails
    assert tan(pi/4)    == 1
#   assert tan(pi/2)    == oo


def _test_asin():
    # XXX No inverse trig
    assert asin(0)      == 0
    assert asin(-1)     == -pi/2

    assert feq( float(asin(Real("0.3"))), 0.30469265401539752 )
    assert asin(Real("0.3")).evalf(precision=28) == Real("0.3046926540153975079720029612")

def _test_acos():
    # XXX No inverse trig
    assert acos(0)      == pi/2
    assert acos(1)      == 0
    assert acos(-1)     == pi

    assert feq( float(acos(Real("0.3"))), 1.2661036727794992 )
    assert acos(Real("0.3")).evalf(precision=28) == Real("1.266103672779499111259318730")

def _test_atan():
    # XXX No inverse trig
    assert atan(1)      == pi/4
    assert atan(oo)     == pi/2
    #assert atan(-oo)    == -pi/2   XXX infinite recursion
    assert atan(-sqrt(3))== - pi/3

    assert feq( float(atan(1)), 0.78539816339744828 )
    assert atan(1).evalf(precision=28) == Real("0.7853981633974483096156608458")
