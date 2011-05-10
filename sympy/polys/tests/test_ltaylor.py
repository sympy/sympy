"""Tests for distributed polynomials and series using lpoly"""
from sympy.polys.ltaylor import taylor

from sympy.functions.elementary.trigonometric import (cos,sin,atan,acos,acot)
from sympy.functions.elementary.exponential import (exp,log,LambertW)
from sympy.functions.elementary.hyperbolic import (asinh,acosh,acoth)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.series.order import O
from sympy.series.series import series
from sympy.core.numbers import (Rational,E)
from sympy.core.symbol import Symbol
from sympy import I
from sympy.core import pi
from sympy.abc import x, y


one = Rational(1)
def test_taylor_QQ1():
    p = x
    p1 = taylor(p,x,0,7)
    assert p1 == x + O(x**7)
    p = -x
    p1 = taylor(p,x,0,7)
    assert p1 == -x + O(x**7)
    p = x+1
    p1 = taylor(p,x,0,11)
    assert p1 == x + 1 + O(x**11)
    p = (1 + x)**20
    p1 = taylor(p,x,0,11)
    assert p1 == 184756*x**10 + 167960*x**9 + 125970*x**8 + 77520*x**7 + 38760*x**6 + 15504*x**5 + 4845*x**4 + 1140*x**3 + 190*x**2 + 20*x + 1 + O(x**11)

    p = 1/(x+1)**2
    p1 = taylor(p,x,0,11)
    assert p1 == 11*x**10 - 10*x**9 + 9*x**8 - 8*x**7 + 7*x**6 - 6*x**5 + 5*x**4 - 4*x**3 + 3*x**2 - 2*x + 1 + O(x**11)

def test_taylor_QQ2():
    p = cos(cos(cos(x) - 1)-1) + cos(cos(x) - 1)
    p1 = taylor(p,x,0,11)
    assert p1 == 437*one/241920*x**10 - 13*one/1920*x**8 + one/48*x**6 - one/8*x**4 + 2 + O(x**11)

    p = 1/cos(x)
    p1 = taylor(p,x,0,11)
    assert p1 == 50521*one/3628800*x**10 + 277*one/8064*x**8 + 61*one/720*x**6 + 5*one/24*x**4 + one/2*x**2 + 1 + O(x**11)

    p = 1/(1+cos(x))
    p1 = taylor(p,x,0,11)
    assert p1 == 691*one/14515200*x**10 + 31*one/80640*x**8 + 17*one/5760*x**6 + one/48*x**4 + one/8*x**2 + one/2 + O(x**11)

    p = cos(cos(x) - 1)/(1+cos(x))
    p1 = taylor(p,x,0,11)
    assert p1 == -283*one/1814400*x**10 + 73*one/80640*x**8 - 13*one/5760*x**6 - one/24*x**4 + one/8*x**2 + one/2 + O(x**11)

def test_taylor_QQ3():
    p = (1 + x)**Rational(1,4)
    p1 = taylor(p,x,0,11)
    assert p1 == -3129357*one/268435456*x**10 + 447051*one/33554432*x**9 - 129789*one/8388608*x**8 + 4807*one/262144*x**7 - 1463*one/65536*x**6 + 231*one/8192*x**5 - 77*one/2048*x**4 + 7*one/128*x**3 - 3*one/32*x**2 + one/4*x + 1 + O(x**11)

    p = (1 + x)**Rational(3,4)
    p1 = taylor(p,x,0,11)
    assert p1 == -1057485*one/268435456*x**10 + 160225*one/33554432*x**9 - 49725*one/8388608*x**8 + 1989*one/262144*x**7 - 663*one/65536*x**6 + 117*one/8192*x**5 - 45*one/2048*x**4 + 5*one/128*x**3 - 3*one/32*x**2 + 3*one/4*x + 1 + O(x**11)

    p = (cos(x))**Rational(-3,4)
    p1 = taylor(p,x,0,11)
    assert p1 == 8787601*one/1238630400*x**10 + 63587*one/3440640*x**8 + 751*one/15360*x**6 + 17*one/128*x**4 + 3*one/8*x**2 + 1 + O(x**11)

def test_taylor_QQ4():
    p = (exp(x))**Rational(-3,4)
    p1 = taylor(p,x,0,7)
    assert p1 == 81*one/327680*x**6 - 81*one/40960*x**5 + 27*one/2048*x**4 - 9*one/128*x**3 + 9*one/32*x**2 - 3*one/4*x + 1 + O(x**7)

    p = (1 + log(1+x))**Rational(-3,4)
    p1 = taylor(p,x,0,7)
    assert p1 == 5114923*one/983040*x**6 - 139809*one/40960*x**5 + 4619*one/2048*x**4 - 193*one/128*x**3 + 33*one/32*x**2 - 3*one/4*x + 1 + O(x**7)

    p = sqrt((1 + log(1+x))**Rational(-3,4))
    p1 = taylor(p,x,0,7)
    assert p1 == 103245539*one/62914560*x**6 - 1497009*one/1310720*x**5 + 26491*one/32768*x**4 - 601*one/1024*x**3 + 57*one/128*x**2 - 3*one/8*x + 1 + O(x**7)

    p = asinh(x**2)
    p1 = taylor(p,x,0,11)
    assert p1 == 3*one/40*x**10 - one/6*x**6 + x**2 + O(x**11)

def test_taylor_QQ_par1():
    p = (x+y)**4
    p1 = taylor(p,x,0,3,pol_pars=[y])
    assert p1 == 6*x**2*y**2 + 4*x*y**3 + y**4 + O(x**3)
    p = (1+x+y)**4
    p1 = taylor(p,x,0,3,pol_pars=[y])
    assert p1 == 6*x**2*y**2 + 4*x*y**3 + y**4 + 12*x**2*y + 12*x*y**2 + 6*x**2 + 4*y**3 + 12*x*y + 6*y**2 + 4*x + 4*y + 1 + O(x**3)

    p = atan(x*y + x**2)
    p1 = taylor(p,x,0,6)
    assert p1.expand() == one/5*x**5*y**5 - x**5*y - x**4*y**2 - one/3*x**3*y**3 + x**2 + x*y + O(x**6)

def test_taylor_SR1():
    x = Symbol('x')
    p = cos(sqrt(2)*x)
    p1 = taylor(p,x,0,11)
    assert p1 == -one/113400*x**10 + one/2520*x**8 - one/90*x**6 + one/6*x**4 - x**2 + 1 + O(x**11)

    p = sqrt(1 + cos(x))
    p1 = taylor(p,x,0,11)
    assert p1 == -one/3715891200*sqrt(2)*x**10 + one/10321920*sqrt(2)*x**8 - one/46080*sqrt(2)*x**6 + one/384*sqrt(2)*x**4 - one/8*sqrt(2)*x**2 + sqrt(2) + O(x**11)

    p = cos(x+40)
    p1 = taylor(p,x,0,5)
    assert p1 == -x*sin(40) - x**2*cos(40)/2 + x**3*sin(40)/6 + x**4*cos(40)/24 + cos(40) + O(x**5)

    p = 1/cos(x+40)
    p1 = taylor(p,x,0,3)
    p2 = x**2*sin(40)**2/cos(40)**3 + one/2*x**2/cos(40) + x*sin(40)/cos(40)**2 + 1/cos(40) + O(x**3)
    assert p1.expand() == p2

    p = sin(x+40)
    p1 = taylor(p,x,0,6)
    assert p1 == one/120*x**5*cos(40) + one/24*x**4*sin(40) - one/6*x**3*cos(40) - one/2*x**2*sin(40) + x*cos(40) + sin(40) + O(x**6)


def test_taylor_SR2():
    p = exp(2 - exp(1+x+x**2))
    p1 = taylor(p,x,0,5)
    p2 = -25*one/24*x**4*E**(-E + 3) + 55*one/24*x**4*E**(-E + 4) - 3*one/4*x**4*E**(-E + 5) + one/24*x**4*E**(-E + 6) - 7*one/6*x**3*E**(-E + 3) + 3*one/2*x**3*E**(-E + 4) - one/6*x**3*E**(-E + 5) - 3*one/2*x**2*E**(-E + 3) + one/2*x**2*E**(-E + 4) - x*E**(-E + 3) + E**(-E + 2) + O(x**5)
    assert p2.expand() == p1.expand()

def test_taylor_series():
    x = Symbol('x')
    p = exp(x*log(x))
    #f = open('/home/pernici/AA/Python/home/rmpoly/lpoly','w')
    p1 = taylor(p,x,0,3)
    assert p1 == 1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)

    p=taylor(1/cos(x*sqrt(x)),x,0,10)
    assert p == 1 + x**3/2 + 5*x**6/24 + 61*x**9/720 + O(x**10)


def test_lambert():
    p = LambertW(x)
    h = 20
    p1 = taylor(p,x,0,h)
    p2 = 0
    fact = 1
    for i in range(1,h):
        fact *= i
        p2 += (-Rational(i))**(i-1)/fact*x**i
    assert p1 == p2 + O(x**h)

def test_acos():
    p = acos(x)
    p1 = taylor(p,x,0,10)
    assert p1 == pi/2 - 35*one/1152*x**9 - 5*one/112*x**7 - 3*one/40*x**5 - one/6*x**3 - x + O(x**10)
    p = acos(x+x**2)
    p1 = taylor(p,x,0,10)
    assert p1 == pi/2 - 1223*one/640*x**10 - 1547*one/1152*x**9 - 17*one/16*x**8 - 89*one/112*x**7 - 13*one/24*x**6 - 23*one/40*x**5 - one/2*x**4 - one/6*x**3 - x**2 - x + O(x**10)

def test_acot():
    p = acot(x+x**2)
    p1 = taylor(p,x,0,10)
    p2 = pi/2 - x - x**2 + x**3/3 + x**4 + 4*x**5/5 - 2*x**6/3 - 13*x**7/7 - x**8 + 17*x**9/9 + O(x**10)
    assert p1 == p2

def test_acoth():
    p = acoth(x+x**2)
    p1 = taylor(p,x,0,10)
    p2 = -I*pi/2 + 37*one/9*x**9 + 3*x**8 + 15*one/7*x**7 + 4*one/3*x**6 + 6*one/5*x**5 + x**4 + one/3*x**3 + x**2 + x + O(x**10)
    assert p1 == p2

def test_acosh():
    p = acosh(x)
    p1 = taylor(p,x,0,10)
    p2 = -I*pi/2 + 35*one/1152*I*x**9 + 5*one/112*I*x**7 + 3*one/40*I*x**5 + one/6*I*x**3 + I*x + O(x**10)
    assert p1 == p2

def test_sum():
    p = log(2) + x
    p1 = taylor(p,x,0,10)
    assert p1 == log(2) + x + O(x**10)
    p = log(2) + sin(x)
    p1 = taylor(p,x,0,8)
    assert p1 == -one/5040*x**7 + one/120*x**5 - one/6*x**3 + x + log(2) + O(x**8)
    p = 1/x + sin(x)
    p1 = taylor(p,x,0,8)
    assert p1 == -one/5040*x**7 + one/120*x**5 - one/6*x**3 + x + one/x + O(x**8)
    p = log(2) + sqrt(3)*sin(x)
    p1 = taylor(p,x,0,8)
    assert p1.expand() == -one/5040*sqrt(3)*x**7 + one/120*sqrt(3)*x**5 - one/6*sqrt(3)*x**3 + sqrt(3)*x + log(2) + O(x**8)

    p = log(2) + (1+sqrt(3))*sin(x) + 1/x
    p1 = taylor(p,x,0,8)
    p2 = -one/5040*(sqrt(3) + 1)*x**7 + one/120*(sqrt(3) + 1)*x**5 - one/6*(sqrt(3) + 1)*x**3 + (sqrt(3) + 1)*x + 1/x + log(2) + O(x**8)
    assert p1.expand() == p2.expand()

