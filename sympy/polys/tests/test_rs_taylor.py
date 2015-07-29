# tests for branch lt1, recovered from branch lpoly2
import sys
from time import time
sys.path.insert(0, '/home/pernici/AA/Python/packages/sympy/branches/fork/sympy')

from sympy import simplify, sqrt, E, Function, cacheit, erf, symbols, Symbol, \
sqrt, pi, erf, expand, atan2, cot
from sympy.polys.domains import QQ, EX
from sympy.polys.rings import ring, xring
from sympy.polys.rs_taylor import *


def test_taylor():
    x, y = symbols('x, y')
    r = taylor(sin(x*tan(x)), x, 0, 10)
    assert r == x**2 + x**4/3 - x**6/30 - 71*x**8/630 + O(x**10)
    p = x*exp(x)/(sin(x)*tan(2*x))
    r = taylor(p, x, 0, 5)
    assert r == 1/(2*x) + S(1)/2 - x/3 - x**2/2 - 11*x**3/20 - 67*x**4/180 + O(x**5)

    r = taylor(atan(x*y + x**2), x, 0, 5)
    assert r == x*y + x**2 - x**3*y**3/3 - x**4*y**2 + O(x**5)
    r = taylor(atan(x*y + x**2), x, 0, 5, pol_pars=[y])
    assert r == x*y + x**2 - x**3*y**3/3 - x**4*y**2 + O(x**5)
    
def test_polynomial_degree():
    x = Symbol('x')
    assert polynomial_degree(S.One, x) == 0
    assert polynomial_degree(x, x) == 1
    assert polynomial_degree(-x, x) == 1
    assert polynomial_degree((x+1)*((x**2+1)**2 + 1)**2 + x + 1, x) == 9
    assert polynomial_degree((x+1)*(x**2+1)**2 - x**5, x) == 5
    assert polynomial_degree(1/x, x) == -1
    assert polynomial_degree(x/(1+x), x) == -1
    assert polynomial_degree(sin(x), x) == -1

def test_poly_truncate():
    x = Symbol('x')
    assert poly_truncate(x.as_poly(x), x, 6) == x
    assert poly_truncate(((1 + x)**4).as_poly(), x, 3) == 6*x**2 + 4*x + 1

def test_taylor_QQ1():
    x = Symbol('x')
    assert taylor(x, x, 0, 7) == x
    assert taylor(-x, x, 0, 7) == -x
    assert taylor(x+1, x, 0, 11) == x+1
    assert taylor(x**4, x, 0, 4) == O(x**4)
    assert taylor(x**4, x, 0, 5) == x**4
    assert taylor((x**2 + 1)**2 - x**4, x, 0, 4) == 2*x**2 + 1
    assert taylor(sqrt(x) + 1, x, 0, 5) == sqrt(x) + 1
    assert taylor(x + O(x**6), x, 0, 4) == x + O(x**4)
    assert taylor(2 + x**10, x, 0, 4) == 2 + O(x**4)
    assert taylor(log(x)*x, x, 0, 1) == x*log(x)
    assert taylor(log(x)*x**3, x, 0, 4) == x**3*log(x)
    assert taylor(log(x)*x**4, x, 0, 4) == x**4*log(x)
    assert taylor(log(x)*x**6, x, 0, 4) == O(x**4)
    assert taylor((x**4 + x**2)*log(x), x, 0, 4) == x**4*log(x) + x**2*log(x)
    assert taylor((x**2 + 1)**2*log(x), x, 0, 4) == x**4*log(x) + 2*x**2*log(x) + log(x)
    assert taylor(x**2*(x**3 + x**2 + 1)*log(x), x, 0, 4) == x**2*log(x) + x**4*log(x) + O(x**4)
    assert taylor((1 + x + x**2 + x**4)**2*log(x), x, 0, 4) == \
      log(x) + 2*x*log(x) + 3*x**2*log(x) + 2*x**3*log(x) + 3*x**4*log(x) + O(x**4)
    assert taylor((1 + x)**20, x, 0, 11) == \
      1 + 20*x + 190*x**2 + 1140*x**3 + 4845*x**4 + 15504*x**5 + 38760*x**6 + 77520*x**7 + 125970*x**8 + 167960*x**9 + 184756*x**10 + O(x**11)

    assert taylor(1/(x+1)**2, x, 0, 11) == \
      1 - 2*x + 3*x**2 - 4*x**3 + 5*x**4 - 6*x**5 + 7*x**6 - 8*x**7 + 9*x**8 - 10*x**9 + 11*x**10 + O(x**11)
    y = Symbol('y', positive=True)
    assert taylor(x**y, x, 0, 4) == x**y
    assert taylor(x**-y, x, 0, 4) == exp(-y*log(x))

    r = taylor((S(8)/27 + sin(x))**Rational(1,3), x, 0, 3)
    assert r == \
      Rational(2,3) + 3*x/4 - 27*x**2/32 + O(x**3)
    r = taylor((1 + 3*x + 2*x**2)**100, x, 0, 5)

    assert taylor((1 + 3*x + 2*x**2)**100, x, 0, 5) == \
        1 + 300*x + 44750*x**2 + 4425300*x**3 + 326370825*x**4 + O(x**5)

def test_taylor_QQ2():
    x = Symbol('x')
    assert taylor(cos(cos(cos(x) - 1)-1) + cos(cos(x) - 1), x, 0, 11) == \
       2 - x**4/8 + x**6/48 - 13*x**8/1920 + 437*x**10/241920 + O(x**11)

    assert taylor(1/cos(x), x, 0, 11) == \
      1 + x**2/2 + 5*x**4/24 + 61*x**6/720 + 277*x**8/8064 + 50521*x**10/3628800 + O(x**11)

    assert taylor(1/(1+cos(x)), x, 0, 11) == \
      S.Half + x**2/8 + x**4/48 + 17*x**6/5760 + 31*x**8/80640 + 691*x**10/14515200 + O(x**11)

    assert taylor(cos(cos(x) - 1)/(1+cos(x)), x, 0, 11) == \
      S.Half + x**2/8 - x**4/24 - 13*x**6/5760 + 73*x**8/80640 - 283*x**10/1814400 + O(x**11)

    assert taylor(1/(sin(x)*cos(x)), x, 0, 5) == 1/x + 2*x/3 + 14*x**3/45 + O(x**5)

    assert taylor(atanh(x)*sin(x), x, 0, 7) == x**2 + x**4/6 + 11*x**6/72 + O(x**7)

def test_taylor_QQ3():
    x = Symbol('x')
    assert taylor((1 + x)**Rational(1,4), x, 0, 11) == \
      1 + x/4 - 3*x**2/32 + 7*x**3/128 - 77*x**4/2048 + 231*x**5/8192 - 1463*x**6/65536 + 4807*x**7/262144 - 129789*x**8/8388608 + 447051*x**9/33554432 - 3129357*x**10/268435456 + O(x**11)

    assert taylor((1 + x)**Rational(3,4), x, 0, 11) == \
      1 + 3*x/4 - 3*x**2/32 + 5*x**3/128 - 45*x**4/2048 + 117*x**5/8192 - 663*x**6/65536 + 1989*x**7/262144 - 49725*x**8/8388608 + 160225*x**9/33554432 - 1057485*x**10/268435456 + O(x**11)

    assert taylor((cos(x))**Rational(-3,4), x, 0, 11) == \
      1 + 3*x**2/8 + 17*x**4/128 + 751*x**6/15360 + 63587*x**8/3440640 + 8787601*x**10/1238630400 + O(x**11)

def test_taylor_QQ4():
    x = Symbol('x')
    assert taylor((exp(x))**Rational(-3,4), x, 0, 7) == \
      1 - 3*x/4 + 9*x**2/32 - 9*x**3/128 + 27*x**4/2048 - 81*x**5/40960 + 81*x**6/327680 + O(x**7)

    assert taylor((1 + log(1+x))**Rational(-3,4), x, 0, 7) == \
      1 - 3*x/4 + 33*x**2/32 - 193*x**3/128 + 4619*x**4/2048 - 139809*x**5/40960 + 5114923*x**6/983040 + O(x**7)

    assert taylor(sqrt((1 + log(1+x))**Rational(-3,4)), x, 0, 7) == \
      1 - 3*x/8 + 57*x**2/128 - 601*x**3/1024 + 26491*x**4/32768 - 1497009*x**5/1310720 + 103245539*x**6/62914560 + O(x**7)

    assert taylor(asinh(x**2), x, 0, 11) == x**2 - x**6/6 + 3*x**10/40 + O(x**11)
    assert taylor(asin(x + x**2), x, 0, 6) == x + x**2 + x**3/6 + x**4/2 + 23*x**5/40 + O(x**6)


def test_taylor_QQ5():
    x = Symbol('x')
    assert taylor(1/(exp(x)-1), x, 0, 7) == 1/x - S.Half + x/12 - x**3/720 + x**5/30240 + O(x**7)

    assert taylor(sin(x)/sin(tan(x + x**2)), x, 0, 7) == \
      1 - x + 2*x**2/3 - x**3 + 83*x**4/90 - 59*x**5/90 + 1949*x**6/1890 + O(x**7)

def test_taylor_QQ6():
    x = Symbol('x')
    assert taylor((1 + 1/x)/cos(x), x, 0, 7) == \
      1/x + 1 + x/2 + x**2/2 + 5*x**3/24 + 5*x**4/24 + 61*x**5/720 + 61*x**6/720 + O(x**7)

    #assert taylor(cos(x)**sin(x), x, 0, 10) == \
    #  1 - x**3/2 + x**6/8 - x**7/80 - 37*x**9/1512 + O(x**10)

    assert taylor((2/x+3/x**2)/(1/x+1/x**2), x, 0, 5) == \
      3 - x + x**2 - x**3 + x**4 + O(x**5)

    assert taylor((exp(x)-1)/x, x, 0, 7) == \
      1 + x/2 + x**2/6 + x**3/24 + x**4/120 + x**5/720 + x**6/5040 + O(x**7)

    assert taylor((exp(x)-1)/(x + 1/x), x, 0, 7) == \
      x**2 + x**3/2 - 5*x**4/6 - 11*x**5/24 + 101*x**6/120 + O(x**7)

    assert taylor((exp(x)-1)/(x + sin(x)/x), x, 0, 7) == \
      x - x**2/2 + 5*x**3/6 - 7*x**4/8 + 73*x**5/72 - 277*x**6/240 + O(x**7)

    assert taylor((exp(x)-1)/(-1 + sin(x)/x), x, 0, 5) == \
      -6/x - 3 - 13*x/10 - 2*x**2/5 - 151*x**3/1400 - 13*x**4/525 + O(x**5)

    assert taylor(sqrt(x*sin(x)), x, 0, 7) == x - x**3/12 + x**5/1440 + O(x**7)
    assert taylor(sqrt(atanh(x)*sin(x)), x, 0, 7) == x + x**3/12 + 7*x**5/96 + O(x**7)

    assert taylor(Pow(x*atanh(x)*sin(x), S.One/3), x, 0, 7) == \
      x + x**3/18 + 31*x**5/648 + O(x**7)


def test_taylor_QQ_logx():
    x = Symbol('x')
    assert taylor(sin(x)**2*log(x), x, 0, 7) == \
      x**2*log(x) - x**4*log(x)/3 + 2*x**6*log(x)/45 + O(x**7)
    # series gives O(x**7*log(x))

    p1 = taylor(sin(x)**2*log(sin(x) + x)*log(x), x, 0, 7)
    p2 = x**2*(log(x)**2 + log(2)*log(x)) + x**4*(-log(x)**2/3 - log(2)*log(x)/3 - log(x)/12) + x**6*(2*log(x)**2/45 + 41*log(x)/1440 + 2*log(2)*log(x)/45) + O(x**7)
    assert p1 == p2
    # series gives O(x**7*log(x)**2)

    p = log(x)*log(x + sin(x))*log(x**2 + sin(2*x))*sin(x)**2*log(x + x**2)
    p1 = taylor(p, x, 0, 4)
    p2 = x**2*(log(x)**4 + 2*log(2)*log(x)**3 + log(2)**2*log(x)**2) + x**3*(3*log(x)**3/2 + 5*log(2)*log(x)**2/2 + log(2)**2*log(x)) + x**4*(-log(x)**4/3 - 11*log(x)**3/8 - 2*log(2)*log(x)**3/3 - 15*log(2)*log(x)**2/8 - log(2)**2*log(x)**2/3 + log(x)**2/2 - log(2)**2*log(x)/2 + log(2)*log(x)/2) + O(x**4)

    assert p1 == p2

    p1 = taylor(tan(x)**2*log(tan(x)), x, 0, 8)
    p2 = x**2*log(x) + x**4*(2*log(x)/3 + S.One/3) + x**6*(17*log(x)/45 + S(3)/10) + 62*x**8*log(x)/315 + O(x**8)

    # the series method gives  order O(x**8*log(x))
    assert p1 == p2

    assert taylor(log(sin(2*x)), x, 0, 7) == \
      log(2) + log(x) - 2*x**2/3 - 4*x**4/45 - 64*x**6/2835 + O(x**7)

    p1 = taylor(log(2/x)*log(sin(2*x)), x, 0, 3)
    p2 = log(2)**2 - log(x)**2 + x**2*(2*log(x)/3 - 2*log(2)/3) + O(x**3)

    #  the series method gives  order O(x**3*log(x))
    assert p1 == p2

    p1 = taylor(log(2/x)*log(sin(2*x)), x, 0, 4)
    assert p1 == log(2)**2 - log(x)**2 + x**2*(2*log(x)/3 - 2*log(2)/3) + 4*x**4*log(x)/45 + O(x**4)

    # series has order O(x**4*log(x)), missing last term
    y = Symbol('y', positive=True)
    assert taylor(x**y*log(x), x, 0, 4) == x**y*log(x)

    p = taylor(1/cos(x*sqrt(x)), x, 0, 10)
    assert p == 1 + x**3/2 + 5*x**6/24 + 61*x**9/720 + O(x**10)

    p1 = taylor(cos(sqrt(2)*x)*log(tan(x)), x, 0, 6)
    p2 = log(x) + x**2*(-log(x) + S.One/3) + x**4*(log(x)/6 - S(23)/90) - x**6*log(x)/90 + O(x**6)
    assert p1 == p2

def test_taylor_QQ_root():
    x = Symbol('x')
    assert taylor(1/Pow(x*sin(x), S.Half), x, 0, 7) == \
      1/x + x/12 + x**3/160 + 61*x**5/120960 + O(x**7)

    assert taylor(Pow(x*sin(x), S.One/2), x, 0, 9) == \
      x - x**3/12 + x**5/1440 - x**7/24192 + O(x**9)

    assert taylor(Pow(x*sin(x)*tan(x), S.One/3), x, 0, 9) == \
      x + x**3/18 + 83*x**5/3240 + 10457*x**7/1224720 + O(x**9)

    assert taylor(Pow(x*sin(x)*tan(x), -S.One/3), x, 0, 7) == \
      1/x - x/18 - 73*x**3/3240 - 7181*x**5/1224720 + O(x**7)

    assert taylor(sin(x)/(cos(x)**2 * sqrt(x*sin(x))), x, 0, 7) == \
    1 + 11*x**2/12 + 841*x**4/1440 + 7811*x**6/24192 + O(x**7)

    p1 = taylor(sin(x)/(cos(x)**2 * Pow(x*sin(x)*tan(x), S.One/3)), x, 0, 8)
    assert p1 == 1 + 7*x**2/9 + 178*x**4/405 + 16987*x**6/76545 + O(x**8)

    assert taylor(sin(x)*Pow(x*sin(x)*tan(x), S.One/3)/cos(x)**2, x, 0, 9) == \
      x**2 + 8*x**4/9 + 47*x**6/81 + 25484*x**8/76545 + O(x**9)

    p1 = taylor(sin(x)*Pow(x*sin(x)*tan(x), S.One/3)/cos(x)**2, x, 0, 9)
    assert p1 == x**2 + 8*x**4/9 + 47*x**6/81 + 25484*x**8/76545 + O(x**9)

    p1 = taylor(sin(x)*Pow((1+cos(x))*x*sin(x)*tan(x), S.One/3)/cos(x)**2, x, 0, 8)
    c = Pow(2, S.One/3)
    assert p1 == c*x**2 + 29*c*x**4/36 + 41*c*x**6/81 + O(x**8)

    p1 = taylor(sin(x)/(cos(x)**2 * Pow((1+cos(x))*x*sin(x)*tan(x), S.One/3)), x, 0, 5)
    c = Pow(2, S(2)/3)
    assert p1 == c/2 + 31*c*x**2/72 + 3313*c*x**4/12960 + O(x**5)

def test_taylor_SR1():
    x = Symbol('x')
    p1 = taylor(cos(sqrt(2)*x), x, 0, 11)
    assert p1 == 1 - x**2 + x**4/6 - x**6/90 + x**8/2520 - x**10/113400 + O(x**11)

    p1 = taylor(sqrt(1 + cos(x)), x, 0, 11)
    p2 = sqrt(2)*(1 - x**2/8 + x**4/384 - x**6/46080 + x**8/10321920 - x**10/3715891200 + O(x**12)) + O(x**11)
    # ATT series method collects sqrt(2)
    assert p1 == p2.expand()

    p1 = taylor(cos(x+40), x, 0, 5)
    assert p1 == cos(40) - x*sin(40) - x**2*cos(40)/2 + x**3*sin(40)/6 + x**4*cos(40)/24 + O(x**5)

    p1 = taylor(1/cos(x+40), x, 0, 3)
    p2 = 1/cos(40) + x*sin(40)/cos(40)**2 + x**2*(sin(40)**2/cos(40)**3 + 1/(2*cos(40))) + O(x**3)
    p1 = expand(p1)
    p2 = expand(p2)
    assert p1 == p2

    """
    p1 = taylor(sin(x+40), x, 0, 6)
    # FIXME
    assert p1 == sin(40) + x*cos(40) - x**2*sin(40)/2 - x**3*cos(40)/6 + x**4*sin(40)/24 + x**5*cos(40)/120 + O(x**6)
    """

    p1 = taylor(sin(x)/((1 + cos(x)/2)*(2 + cos(x)/2)**S.Half), x, 0, 4)
    assert p1 == 2*sqrt(10)*x/15 + sqrt(10)*x**3/150 + O(x**4)

    p1 = taylor(1/(cos(x)*(cos(2) + x)), x, 0, 2)
    p1 = expand(p1)
    assert p1 == 1/cos(2) - x/cos(2)**2 + O(x**2)

    p1 = taylor(1/((cos(2) + x)*(1 + cos(x))), x, 0, 2)
    p1 = expand(p1)
    assert p1 == 1/(2*cos(2)) - x/(2*cos(2)**2) + O(x**2)

    p1 = taylor(1/(cos(x/3)*(cos(2) + x/7)), x, 0, 2)
    p1 = expand(p1)
    assert p1 == 1/cos(2) - x/(7*cos(2)**2) + O(x**2)

def test_taylor_SR2():
    x = Symbol('x')
    p1 = taylor(exp(2 - exp(1+x+x**2)), x, 0, 5)
    p2 = exp(-E + 2) - x*exp(-E + 3) + x**2*exp(-E + 4)/2 - 3*x**2*exp(-E + 3)/2 + 3*x**3*exp(-E + 4)/2 - 7*x**3*exp(-E + 3)/6 - x**3*exp(-E + 5)/6 + 55*x**4*exp(-E + 4)/24 + x**4*exp(-E + 6)/24 - 25*x**4*exp(-E + 3)/24 - 3*x**4*exp(-E + 5)/4 + O(x**5)
    # ATT series collects exponentials
    assert p1.expand() == p2.expand()

    assert taylor(sin((1 + I) + x), x, 0, 2) == sin(1 + I) + x*cos(1 + I) + O(x**2)

def test_SR3():
    # example from http://groups.google.com/group/sympy/msg/37fea5ab11302a08
    x,y = symbols('x,y')
    p = (exp(x))/((y/2 + log(2*pi)/2 + x/12 - 1/x - y/x))
    p1 = taylor(taylor(p, x, 0, 3).removeO(), y, 0, 2)
    p2 = -x + x*y - x**2*log(2)/2 - x**2*log(pi)/2 + y*x**2*log(2) + y*x**2*log(pi) - x**2 + y*x**2/2 + O(y**2)

    assert p1 == p2

    p1 = taylor(log(1 - cos(x**4)), x, 0, 1)
    assert p1 == -log(2) + 8*log(x) + O(x)

def test_taylor_SR_logx():
    x = Symbol('x')
    p1 = taylor(sin(x)**2*log(sin(x) + pi*x)*log(x), x, 0, 6)
    p2 = x**2*log(x)*log(1 + pi) + x**2*log(x)**2 - 120*x**4*log(x)/(720 + 720*pi) - x**4*log(x)*log(1 + pi)/3 - x**4*log(x)**2/3 + O(x**6*log(x)**2)
    assert simplify(p1 - p2) == O(x**6*log(x)**2)

    p1 = taylor(tan(pi*x)**2*log(tan(x)), x, 0, 6)
    p2 = pi**2*x**2*log(x) + x**4*(2*pi**4*log(x)/3 + pi**2/3) + 17*pi**6*x**6*log(x)/45 + O(x**6)
    # series has order O(x**6*log(x)), missing last term
    assert p1 == p2

def test_taylor_SR_root():
    x = Symbol('x')
    p1 = taylor(sqrt(sin(cos(2)*x)), x, 0, 7)
    assert p1 == sqrt(x)*sqrt(cos(2)) - x**(S(5)/2)*cos(2)**(S(5)/2)/12 + x**(S(9)/2)*cos(2)**(S(9)/2)/1440 - x**(S(13)/2)*cos(2)**(S(13)/2)/24192 + O(x**7)

    p1 = taylor(sqrt(x*sin(cos(2)*x)), x, 0, 7)
    c = cos(2)
    assert p1 == x*sqrt(c) - x**3*c**(S(5)/2)/12 + x**5*c**(S(9)/2)/1440 + O(x**7)

def notest_lambert():
    x = Symbol('x')
    p = LambertW(x)
    h = 20
    p1 = taylor(LambertW(x), x, 0, h)
    p2 = 0
    fact = 1
    for i in range(1,h):
        fact *= i
        p2 += (-Rational(i))**(i-1)/fact*x**i
    assert p1 == p2 + O(x**h)

def test_acos():
    x = Symbol('x')
    p1 = taylor(acos(x), x, 0, 10)
    assert p1 == pi/2 - x - x**3/6 - 3*x**5/40 - 5*x**7/112 - 35*x**9/1152 + O(x**10)

    p1 = taylor(acos(x+x**2), x, 0, 10)
    assert p1 == pi/2 - x - x**2 - x**3/6 - x**4/2 - 23*x**5/40 - 13*x**6/24 - 89*x**7/112 - 17*x**8/16 - 1547*x**9/1152 + O(x**10)

def test_acot():
    x = Symbol('x')
    p1 = taylor(acot(x+x**2), x, 0, 10)
    p2 = pi/2 - x - x**2 + x**3/3 + x**4 + 4*x**5/5 - 2*x**6/3 - 13*x**7/7 - x**8 + 17*x**9/9 + O(x**10)
    assert p1 == p2

def test_acoth():
    x = Symbol('x')
    p1 = taylor(acoth(x+x**2), x, 0, 10)
    p2 = I*pi/2 + x + x**2 + x**3/3 + x**4 + 6*x**5/5 + 4*x**6/3 + 15*x**7/7 + 3*x**8 + 37*x**9/9 + O(x**10)
    assert p1 == p2

def test_acosh():
    x = Symbol('x')
    p1 = taylor(acosh(x), x, 0, 10)
    p2 = I*pi/2 - I*x - I*x**3/6 - 3*I*x**5/40 - 5*I*x**7/112 - 35*I*x**9/1152 + O(x**10)
    assert p1 == p2

def test_sum():
    x = Symbol('x')
    p1 = taylor(log(2) + x, x, 0, 10)
    assert p1 == log(2) + x

    p1 = taylor(log(2) + sin(x), x, 0, 8)
    assert p1 == log(2) + x - x**3/6 + x**5/120 - x**7/5040 + O(x**8)

    p = 1/x + sin(x)
    p1 = taylor(p,x,0,8)
    assert p1 == 1/x + x - x**3/6 + x**5/120 - x**7/5040 + O(x**8)
    p = log(2) + sqrt(3)*sin(x)
    p1 = taylor(p,x,0,8)
    assert p1 == log(2) + sqrt(3)*x - sqrt(3)*x**3/6 + sqrt(3)*x**5/120 - sqrt(3)*x**7/5040 + O(x**8)

    p1 = taylor(log(2) + (1+sqrt(3))*sin(x) + 1/x, x, 0, 8)
    p2 = 1/x + log(2) + sqrt(3)*x + x - x**3/6 - sqrt(3)*x**3/6 + sqrt(3)*x**5/120 + x**5/120 - x**7/5040 - sqrt(3)*x**7/5040 + O(x**8)
    assert p1 == p2

def test_factor_var_from_num():
    x = Symbol("x")
    assert taylor(1 + 1/x**2, x, 0, 10) == 1 + 1/x**2

    p = (1 + 1/x**2)*(1 + cos(x)/x)
    p1 = taylor(p,x, 0,5)
    assert p1 == 1 - 11*x/24 + 1/(2*x) + x**(-3) + x**(-2) + 29*x**3/720 + O(x**5)
    p = (1 + 1/x**2)*(1 + sqrt(2)*cos(x)/x)
    p1 = taylor(p,x, 0,5)
    sq2 = 2**Rational(1,2)
    assert p1 == 1 - 11*x*sq2/24 + sq2/(2*x) + sq2/x**3 + x**(-2) + 29*sq2*x**3/720 + O(x**5)

    p = (1 + 1/x**2)*(1 + sqrt(2)*cos(x)/x)/sin(x)
    p1 = taylor(p,x, 0,5)
    assert p1 == 67*x/360 - 16*sq2/45 + 7/(6*x) + sq2/x**4 + x**(-3) + 2*sq2/(3*x**2) - 23*sq2*x**2/945 + 65*x**3/3024 - 11*sq2*x**4/4725 + O(x**5)

def test_factor_var_from_num_1():
    x = Symbol("x")
    p = ((1 + 1/x)*(1 + 2/x))
    p1 = taylor(((1 + 1/x)*(1 + 2/x)), x, 0, 5)
    # the series method does not give the order, taylor does
    assert p1 == 1 + 3/x + 2/x**2 + O(x**5)

def test_singular():
    x = Symbol('x')
    p1 = taylor(1/(1 + sin(x)/x**2), x, 0, 5)
    assert p1 == x - x**2 + 7*x**3/6 - 4*x**4/3 + O(x**5)

    p1 = taylor(1/(1 + sin(x)/x**2)**2, x, 0, 5)
    assert p1 == x**2 - 2*x**3 + 10*x**4/3 + O(x**5)

    p1 = taylor(1/(1 + sin(x)/x**2)**4, x, 0, 5)
    assert p1 == x**4 + O(x**5)

    p1 = taylor(1/(1 + 1/x + sin(x)/x**2), x, 0, 5)
    assert p1 == x/2 - x**2/4 + x**3/6 - 5*x**4/48 + O(x**5)

    assert taylor((x+x**2)**2, x, 0, 4) == x**2 + 2*x**3 + O(x**4)

    assert taylor((1/x + x**2)**2, x, 0, 4) == x**(-2) + 2*x + O(x**4)

    assert taylor((log(2)+1/x)**2,x,0,4) == log(2)**2+2*log(2)/x+x**(-2)+O(x**4)

    p1 = taylor((1 + (sin(x)/x**2))**3, x, 0, 7)
    assert p1 == -47*x/120 + 5/(2*x) + x**(-3) + 3/x**2 + 2*x**2/15 + 173*x**3/15120 - x**4/105 + 311*x**5/604800 + 2*x**6/4725 + O(x**7)

    y = Symbol('y')
    p1 = taylor(1/(1 + y + 1/x), x, 0, 5)
    assert p1 == x + x**2*(-y - 1) + x**3*(y**2 + 2*y + 1) + x**4*(-y**3 - 3*y**2 - 3*y - 1) + O(x**5)

    p1 = taylor(1/(1 + y/x + 1/x**2), x, 0, 5)
    assert p1 == x**2 - x**3*y + x**4*(y**2 - 1) + O(x**5)

    p1 = taylor(1/(1 + y/x), x, 0, 5)
    assert p1 == x/y - x**4/y**4 + x**3/y**3 - x**2/y**2 + O(x**5)

    p1 = taylor(cos(x)/(x*sin(2*x))**2, x, 0, 7)
    assert p1 == 1/(4*x**4) + 5/(24*x**2) + S(53)/480 + 599*x**2/12096 + 7193*x**4/345600 + 90899*x**6/10644480 + O(x**7)

    p1 = taylor((1 + x**7)/(x + 1/x + exp(x)), x, 0, 10)
    assert p1 == x - x**2 - x**3 + 5*x**4/2 - x**5/6 - 101*x**6/24 + 419*x**7/120 + 4061*x**8/720 - 3557*x**9/336 + O(x**10)

    p1 = taylor(cos(x)/(x + 1/x + exp(sin(x))), x, 0, 10)
    p2 = x - x**2 - 3*x**3/2 + 3*x**4 + 13*x**5/24 - 17*x**6/3 + 2177*x**7/720 + 701*x**8/90 - 433159*x**9/40320 + O(x**10)
    assert p1 == p2

class ftoy1(Function):
    nargs = 1

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
        return x**n

class ftoy2(Function):
    nargs = 1

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n <= 0:
            return S.Zero
        else:
            x = sympify(x)
        return 2*x**-n

class ftoy3(Function):
    nargs = 1

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        x = sympify(x)
        return x**n

def test_analytic():
    x = Symbol('x')
    y = Symbol('y')
    p1 = taylor(log(sin(2*x)*erf(x)*sqrt(pi)), x, 0, 7, analytic=True)
    assert p1 == log(4) + 2*log(x) - x**2 - 2*x**4/45 - 8*x**6/315 + O(x**7)

    p1 = taylor(sqrt(pi)*erf(sin(x)), x, 0, 10, analytic=True)
    assert p1 == 2*x - x**3 + 11*x**5/20 - 241*x**7/840 + 2777*x**9/20160 + O(x**10)
    p1 = taylor(log(sin(2*x)*erf(x)*sqrt(pi)), x, 0, 10, analytic=True)
    assert p1 == log(4) + 2*log(x) - x**2 - 2*x**4/45 - 8*x**6/315 - 4*x**8/567 + O(x**10)

    p1 = taylor(sqrt(pi)*erf(sin(cos(2)*x)), x, 0, 8, analytic=True)
    assert p1 == 2*x*cos(2) - x**3*cos(2)**3 + 11*x**5*cos(2)**5/20 - 241*x**7*cos(2)**7/840 + O(x**8)

    assert taylor(ftoy1(x), x, 0, 1, analytic=True) == O(x)
    assert taylor(ftoy1(x), x, 0, 2, analytic=True) == x + O(x**2)
    assert taylor(ftoy1(x), x, 0, 4, analytic=True) == x + x**3 + O(x**4)

    assert taylor(atan2(x, y), x, 0, 2, analytic=True) == atan2(0, y) + x/y + O(x**2)
    assert taylor(erf(x + 1), x, 0, 2, analytic=True) == erf(1) + 2*x*exp(-1)/sqrt(pi) + O(x**2)
    assert taylor(cot(x), x, 0, 4, analytic=True) == 1/x - x/3 - x**3/45 + O(x**4)
    assert taylor(ftoy2(x), x, 0, 2, analytic=True) == 2/x**3 + 2/x**2 + 2/x + O(x**2)
    assert taylor(ftoy3(x), x, 0, 2, analytic=True) == 1 + x + O(x**2)
    r = taylor(log(sin(2*x)*erf(x)*sqrt(pi)), x, 0, 10, analytic=True)
    assert r == log(4) + 2*log(x) - x**2 - 2*x**4/45 - 8*x**6/315 - 4*x**8/567 + O(x**10)

def test_taylor_series():
    # examples for which currently taylor calls the series method
    x = Symbol('x')
    y = Symbol('y')
    r = taylor(exp(x*log(x)), x, 0, 3)
    assert r ==  1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)
    p1 = taylor(exp(x*log(x)), x, 0, 3)
    assert p1 == 1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)

    assert taylor(2*exp(1/x), x, 0, 4) == 2*exp(1/x)
    # FIXME O(x**2) vs O(x*log(x))
    #assert taylor(log(x)*log(1 - cos(x**4)), x, 0, 1) == \
    #        -log(2)*log(x) + 8*log(x)**2 + O(x*log(x))
    assert taylor((1 + exp(1/x))**2, x, 0, 2) == exp(2/x) + 2*exp(1/x) + 1
    assert taylor((1 + x)**x * x, x, 0, 2) == x + O(x**2)
    assert taylor((1 + 1/sin(x))/sin(x), x, 0, 2) == x**(-2) + 1/x + S.One/3 + x/6 + O(x**2)
    assert taylor(exp(1/x), x, 0, 2) == exp(1/x)
    assert taylor(exp(1/x)*(1 + sin(x)/x), x, 0, 2) == 2*exp(1/x) + O(x**2*exp(1/x))
    assert taylor(sin(x**(1+x)), x, 0, 2) == exp((x + 1)*log(x)) + O(x**2)
