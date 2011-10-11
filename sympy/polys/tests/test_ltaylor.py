"""Tests for distributed polynomials and series using lpoly"""
from sympy.polys.ltaylor import taylor

from sympy.functions.elementary.trigonometric import (cos,sin,atan,acos,acot,tan, asin, cot)
from sympy.functions.elementary.exponential import (exp,log,LambertW)
from sympy.functions.elementary.hyperbolic import (asinh,acosh,acoth,sinh,cosh,tanh)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.series.order import O
from sympy.series.series import series
from sympy.core.numbers import (Rational,E,Integer)
from sympy.core.symbol import Symbol
from sympy import (Symbol, I, Rational,E,Integer,PoleError,floor,ceiling, Piecewise, Eq, sign, symbols, limit, Derivative)
from sympy.core import pi
from sympy.abc import x, y
from sympy.utilities.pytest import raises
from sympy.polys.monomialtools import lex
from sympy import expand, simplify

one = Rational(1)
def test_taylor_QQ1():
    p = x
    p1 = taylor(p,x,0,7)
    assert p1 == x
    p = -x
    p1 = taylor(p,x,0,7)
    assert p1 == -x
    p = x+1
    p1 = taylor(p,x,0,11)
    assert p1 == x + 1
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

    p = series(1/(sin(x)*cos(x)), x, 0,5)
    p1 = taylor(p,x,0,5)
    assert p1 == 1/x + 2*x/3 + 14*x**3/45 + O(x**5)

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

def test_taylor_QQ5():
    p = 1/(exp(x)-1)
    p1 = taylor(p,x,0,7)
    assert p1 == one/30240*x**5 - one/720*x**3 + one/12*x + one/x - one/2 + O(x**7)

    p = sin(x)/sin(tan(x + x**2))
    p1 = taylor(p,x,0,7)
    p2 = 1949*one/1890*x**6 - 59*one/90*x**5 + 83*one/90*x**4 - x**3 + 2*one/3*x**2 - x + 1 + O(x**7)
    assert p1 == p2

def test_taylor_QQ6():
    p = (1 + 1/x)/cos(x)
    p1 = taylor(p,x,0,7)
    assert p1 == 1 + x/2 + 1/x + x**2/2 + 5*x**3/24 + 5*x**4/24 + 61*x**5/720 + 61*x**6/720 + O(x**7)
    p = cos(x)**sin(x)
    p1 = series(p,x, 0,10)
    assert p1 == 1 - x**3/2 + x**6/8 - x**7/80 - 37*x**9/1512 + O(x**10)

    pw=2
    p = (2/x+3/x**pw)/(1/x+1/x**pw)
    p1 = series(p,x,0,5)
    assert p1 == x**4 - x**3 + x**2 - x + 3 + O(x**5)
    p = (exp(x)-1)/x
    p1 = taylor(p,x,0,7)
    assert p1 == 1 + x/2 + x**2/6 + x**3/24 + x**4/120 + x**5/720 + x**6/5040 + O(x**7)
    p = (exp(x)-1)/(x + 1/x)
    p1 = taylor(p,x,0,7)
    assert p1 == x**2 + x**3/2 - 5*x**4/6 - 11*x**5/24 + 101*x**6/120 + O(x**7)
    p = (exp(x)-1)/(x + sin(x)/x)
    p1 = taylor(p,x,0,7)
    #print p.series(x,0,7)
    p2 = x - x**2/2 + 5*x**3/6 - 7*x**4/8 + 73*x**5/72 - 277*x**6/240 + O(x**7)
    assert p1 == p2
    p = (exp(x)-1)/(-1 + sin(x)/x)
    p1 = taylor(p,x,0,5)
    assert p1 == -3 - 13*x/10 - 6/x - 2*x**2/5 - 151*x**3/1400 - 13*x**4/525 + O(x**5)

def test_taylor_QQ_logx():
    p = sin(x)**2*log(x)
    p1 = series(p, x, 0, 7)
    assert p1 == (2*x**6/45 - x**4/3 + x**2)*log(x) + O(x**7*log(x))

    p = sin(x)**2*log(sin(x) + x)*log(x)
    p1 = series(p, x, 0, 7)
    p2 = x**2*log(2)*log(x) + x**2*log(x)**2 - x**4*log(x)/12 - x**4*log(2)*log(x)/3 - x**4*log(x)**2/3 + 2*x**6*log(2)*log(x)/45 + 41*x**6*log(x)/1440 + 2*x**6*log(x)**2/45
    assert expand(p1 - p2) == O(x**7*log(x)**2)

    p = log(x)*log(x + sin(x))*log(x**2 + sin(2*x))*sin(x)**2*log(x + x**2)
    p1 = series(p, x, 0, 4)
    # p2 = p.series(x, 0, 4)
    p2 = x**2*log(2)**2*log(x)**2 + 2*x**2*log(2)*log(x)**3 + x**2*log(x)**4 + x**3*log(2)**2*log(x) + 5*x**3*log(2)*log(x)**2/2 + 3*x**3*log(x)**3/2 + O(x**4*log(x)**4)
    assert expand(p1 - p2) == O(x**4*log(x)**4)

    p = tan(x)**2*log(tan(x))
    p1 = series(p, x, 0, 7)
    # p2 = p.series(x, 0, 7)
    p2 = x**2*log(x) + x**4/3 + 2*x**4*log(x)/3 + 3*x**6/10 + 17*x**6*log(x)/45 + O(x**7*log(x))
    assert expand(p1 - p2) == O(x**7*log(x))

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

def test_SR3():
    # example from http://groups.google.com/group/sympy/msg/37fea5ab11302a08
    x,y = symbols('x,y')
    p = (exp(x))/((y/2 + log(2*pi)/2 + x/12 - 1/x - y/x))
    p1 = series(series(p,x,0,3).removeO(),y,0,2)
    #print 'DBSR3',(p.series(x,0,3).removeO()).series(y,0,2)
    p2 = -x + x*y - x**2*log(2)/2 - x**2*log(pi)/2 + y*x**2*log(2) + y*x**2*log(pi) - x**2 + y*x**2/2 + O(y**2)
    assert p1 == p2

def test_taylor_S_logx():
    p = sin(x)**2*log(sin(x) + pi*x)*log(x)
    p1 = series(p, x, 0, 6)
    # p2 = p.series(x, 0, 6)
    p2 = x**2*log(x)*log(1 + pi) + x**2*log(x)**2 - 120*x**4*log(x)/(720 + 720*pi) - x**4*log(x)*log(1 + pi)/3 - x**4*log(x)**2/3 + O(x**6*log(x)**2)
    assert simplify(p1 - p2) == O(x**6*log(x)**2)

    p = tan(pi*x)**2*log(tan(x))
    p1 = series(p, x, 0, 7)
    # p2 = p.series(x, 0, 7)
    p2 = pi**2*x**2*log(x) + pi**2*x**4/3 + 2*pi**4*x**4*log(x)/3 + 2*pi**4*x**6/9 + 7*pi**2*x**6/90 + 17*pi**6*x**6*log(x)/45 + O(x**7*log(x))
    assert expand(p1 - p2) == O(x**7*log(x))

def test_taylor_series():
    x = Symbol('x')
    p = exp(x*log(x))
    p1 = taylor(p,x,0,3)
    assert p1 == 1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)

    p=taylor(1/cos(x*sqrt(x)),x,0,10)
    assert p == 1 + x**3/2 + 5*x**6/24 + 61*x**9/720 + O(x**10)

    p = cos(sqrt(2)*x)*log(tan(x))
    p1 = series(p, x, 0, 6)
    # p2 = p.series(x, 0, 6)
    p2 = log(x) + x**2/3 - x**2*log(x) - 23*x**4/90 + x**4*log(x)/6 - x**6*log(x)/90 + O(x**6)
    assert expand(p1-p2) == O(x**6*log(x))

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
    # follow series method convention
    p2 = I*pi/2 + 37*one/9*x**9 + 3*x**8 + 15*one/7*x**7 + 4*one/3*x**6 + 6*one/5*x**5 + x**4 + one/3*x**3 + x**2 + x + O(x**10)
    assert p1 == p2

def test_acosh():
    p = acosh(x)
    p1 = taylor(p,x,0,10)
    p2 = -I*pi/2 + 35*one/1152*I*x**9 + 5*one/112*I*x**7 + 3*one/40*I*x**5 + one/6*I*x**3 + I*x + O(x**10)
    # follow series method convention
    assert p1 == -p2

def test_sum():
    p = log(2) + x
    p1 = taylor(p,x,0,10)
    assert p1 == log(2) + x
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


def test_factor_var_from_num():
    x = Symbol("x")
    p = 1 + 1/x**2
    p1 = series(p,x, 0,10)
    assert p1 == 1 + 1/x**2

    p = ((1 + 1/x)*(1 + 2/x))
    #p = cos(x)/x**3
    p1 = series(p,x,0,5)
    # the series method does not give the order
    assert p1 == 1 + 3/x + 2/x**2 + O(x**5)

    p = (1 + 1/x**2)*(1 + cos(x)/x)
    p1 = series(p,x, 0,5)
    assert p1 == 1 - 11*x/24 + 1/(2*x) + x**(-3) + x**(-2) + 29*x**3/720 + O(x**5)
    p = (1 + 1/x**2)*(1 + sqrt(2)*cos(x)/x)
    p1 = series(p,x, 0,5)
    sq2 = 2**Rational(1,2)
    assert p1 == 1 - 11*x*sq2/24 + sq2/(2*x) + sq2/x**3 + x**(-2) + 29*sq2*x**3/720 + O(x**5)

    p = (1 + 1/x**2)*(1 + sqrt(2)*cos(x)/x)/sin(x)
    p1 = series(p,x, 0,5)
    assert p1 == 67*x/360 - 16*sq2/45 + 7/(6*x) + sq2/x**4 + x**(-3) + 2*sq2/(3*x**2) - 23*sq2*x**2/945 + 65*x**3/3024 - 11*sq2*x**4/4725 + O(x**5)

def test_singular():
    x = Symbol('x')
    p = 1/(1 + sin(x)/x**2)
    #p = cos(x)/x**3
    p1 = series(p,x,0,5)
    assert p1 == x - x**2 + 7*x**3/6 - 4*x**4/3 + O(x**5)
    p = 1/(1 + sin(x)/x**2)**2
    p1 = series(p,x,0,5)
    assert p1 == x**2 - 2*x**3 + 10*x**4/3 + O(x**5)
    p = 1/(1 + sin(x)/x**2)**4
    p1 = series(p,x,0,5)
    assert p1 == x**4 + O(x**5)
    p = 1/(1 + 1/x + sin(x)/x**2)
    p1 = series(p,x,0,5)
    assert p1 == x/2 - x**2/4 + x**3/6 - 5*x**4/48 + O(x**5)
    prec = 4
    for p in [(x+x**2)**2,(1/x + x**2)**2,(log(2) + 1/x)**2]:
        p1 = taylor(p,x,0,prec)
        p2 = p.series(x,0,4)
        assert p1 == p2 + O(x**prec)

    p = (1 + (sin(x)/x**2))**3
    p1 = taylor(p,x,0,7)
    assert p1 == -47*x/120 + 5/(2*x) + x**(-3) + 3/x**2 + 2*x**2/15 + 173*x**3/15120 - x**4/105 + 311*x**5/604800 + 2*x**6/4725 + O(x**7)

    y = Symbol('y')
    p = 1/(1 + y + 1/x)
    p1 = series(p,x,0,5)
    assert p1 == x - x**2 + x**3*y**2 - y*x**2 - 3*x**4*y**2 + x**3 - x**4*y**3 + 2*y*x**3 - x**4 - 3*y*x**4 + O(x**5)

    p = 1/(1 + y/x + 1/x**2)
    p1 = series(p,x,0,5)
    assert p1 == x**2 + x**4*y**2 - y*x**3 - x**4 + O(x**5)

    p = 1/(1 + y/x)
    p1 = series(p,x,0,5)
    assert p1 == x/y - x**4/y**4 + x**3/y**3 - x**2/y**2 + O(x**5)

    p = cos(x)/(x*sin(2*x))**2
    p1 = taylor(p,x,0,7)
    assert p1 == 53*one/480 + 1/(4*x**4) + 5/(24*x**2) + 599*x**2/12096 + 7193*x**4/345600 + 90899*x**6/10644480 + O(x**7)

    p = (1 + x**7)/(x + 1/x + exp(x))
    p1 = taylor(p,x,0,10)
    assert p1 == x - x**2 - x**3 + 5*x**4/2 - x**5/6 - 101*x**6/24 + 419*x**7/120 + 4061*x**8/720 - 3557*x**9/336 + O(x**10)

    p = cos(x)/(x + 1/x + exp(sin(x)))
    p1 = taylor(p,x,0,10)
    p2 = x - x**2 - 3*x**3/2 + 3*x**4 + 13*x**5/24 - 17*x**6/3 + 2177*x**7/720 + 701*x**8/90 - 433159*x**9/40320 + O(x**10)
    assert p1 == p2

#######################################################
# tests from test_nseries.py, adapted

def test_nseries():
    """tests adapted from test_nseries.py"""
    assert series(x,x,0,5) == x
    assert series(y,x,0,5) == y

def test_simple_1():
    assert series(x,x,0,5) == x
    y = Symbol("y")
    assert series(y,x,0,5) == y
    p1 = series(1/(x*y),x,0,5)
    assert p1 == 1/(x*y)
    assert series(Rational(3,4),x,0,5) == Rational(3,4)
    assert series(x) == x

def test_mul_0():
    p = x*log(x) + O(x**5)
    assert series(p,x,0,5) == p

def test_mul_1():
    p1 = series(x*log(2+x),x,0,5)
    p2 = x*log(2) + x**2/2 - x**3/8 + x**4/24 + O(x**5)
    assert p1 == p2
    assert series(x*log(1+x),x,0,5) == x**2 - x**3/2 + x**4/3 + O(x**5)

def test_pow_0():
    assert series(x**2,x,0,5) == x**2
    assert series(1/x,x,0,5) == 1/x
    assert series(1/x**2,x,0,5) == 1/x**2
    assert series(x**(Rational(2,3)),x,0,5) == (x**(Rational(2,3)))
    assert series(x**(Rational(3,2)),x,0,5) == (x**(Rational(3,2)))

def test_pow_1():
    assert series((1+x)**2,x,0,5) == 1+2*x+x**2 + O(x**5)

def test_geometric_1():
    assert series(1/(1-x),x,0,5) == 1+x+x**2+x**3+x**4+O(x**5)
    assert series(x/(1-x),x,0,6) == x + x**2 + x**3 + x**4 + x**5 + O(x**6)
    assert series(x**3/(1-x),x,0,8) == x**3 + x**4 + x**5 + x**6 + x**7 + O(x**8)

def test_sqrt_1():
    assert series(sqrt(1+x),x,0,5) == 1+x/2-x**2/8+x**3/16-5*x**4/128+O(x**5)

def test_exp_1():
    assert series(exp(x),x,0,5) == 1+x+x**2/2+x**3/6+x**4/24 + O(x**5)
    assert series(exp(x),x,0,12) == 1+x+x**2/2+x**3/6+x**4/24+x**5/120+  \
               x**6/720+x**7/5040+x**8/40320+x**9/362880+x**10/3628800+  \
               x**11/39916800 + O(x**12)
    assert series(exp(1/x),x,0,5) == exp(1/x)
    assert series(exp(1/(1+x)),x,0,4) == (E*(1-x-13*x**3/6+3*x**2/2)).expand() + O(x**4)
    p1 = series(exp(2+x),x,0,5)
    p2 = (exp(2)*(1+x+x**2/2+x**3/6+x**4/24)).expand() + O(x**5)
    assert series(exp(2+x),x,0,5) == (exp(2)*(1+x+x**2/2+x**3/6+x**4/24)).expand() + O(x**5)

def test_exp_sqrt_1():
    p1 = series(exp(1+sqrt(x)),x,0,3)
    # series calls the series method; p2n is the result from nseries
    #p2n = (exp(1)*(1+sqrt(x)+x/2+sqrt(x)*x/6)).expand() + O(sqrt(x)**3)
    s = sqrt(x)
    p2 = E + E*x/2 + E*s + E*x**2/24 + E*s**3/6 + E*s**5/120 + O(x**3)
    #print exp(1+sqrt(x)).series(x,0,3)
    assert p1 == p2

def test_power_x_x1():
    p1 = series(exp(x*log(x)),x,0,4)
    p2 = 1+x*log(x)+x**2*log(x)**2/2+x**3*log(x)**3/6 + O(x**4*log(x)**4)
    assert p1 == p2

def test_power_x_x2():
    p1 = series((x**x),x,0,4)
    p2 = 1 + x*log(x) + x**2*log(x)**2/2 + x**3*log(x)**3/6 + O(x**4*log(x)**4)
    assert p1 == p2

def test_log_singular1():
    p1 = series(log(1+1/x),x,0,5)
    p2 = x - log(x) - x**2/2 + x**3/3 - x**4/4 + O(x**5)
    assert p1 == p2

def test_log_power1():
    p1 = series(1 / (1/x + x ** (log(3)/log(2))),x,0,5)
    assert p1 == x - x**(2 + log(3)/log(2)) + O(x**5)

#def test_log_series():
#    pass
    # (NO)
    # not done since the series method does not do it
    #l = Symbol('l')
    #p1 = series(1/(1-log(x)),x,0,5)
    #print (1/(1-log(x))).series(x,0,5)
    #assert p1 == 1/(1-l)

def test_log2():
    p1 = series(log(-1/x),x,0,5)
    assert p1 == -log(x) + log(-1)

#def test_log3():
#    pass
    #NO
    # not done since the series method does not do it
    #l = Symbol('l')
    #e = 1/log(-1/x)

def test_series1():
    x = Symbol("x")
    p = sin(x)
    assert series(p,x,0,0) != 0
    assert series(p,x,0,0) == O(1,x)
    assert series(p,x,0,1) == O(x,x)
    assert series(p,x,0,2) == x + O(x**2, x)
    assert series(p,x,0,3) == x + O(x**3, x)
    assert series(p,x,0,4) == x-x**3/6 + O(x**4, x)
    p = (exp(x)-1)/x
    p1 = series(p,x,0,3)
    #p2n = 1+x/2+O(x**2, x)
    assert series(p,x,0,3) == 1 + x/2 + x**2/6 + O(x**3)
    assert series(x,x,0,3) == x

def test_seriesbug1():
    assert series(1/x,x,0,3) == 1/x
    assert series(x + 1/x,x,0,3) == x + 1/x

def test_series2x():
    assert series((x+1)**(-2),x,0,4) == 1-2*x+3*x**2-4*x**3+O(x**4, x)
    assert series((x+1)**(-1),x,0,4) == 1-x+x**2-x**3+O(x**4, x)
    assert series((x+1)**0,x,0,3) == 1
    assert series((x+1)**1,x,0,3) == 1+x
    # O(x**3) should be absent
    assert series((x+1)**2,x,0,3) == 1+2*x+x**2 + O(x**3)
    #p2n = 1 + 3*x + 3*x**2 + x**3
    assert series((x+1)**3,x,0,3) == 1+3*x+3*x**2 + O(x**3)
    assert series(1/(1+x),x,0,4) == 1-x+x**2-x**3+O(x**4, x)
    assert series(x+3/(1+2*x),x,0,4) == 3-5*x+12*x**2-24*x**3+O(x**4, x)
    # O(x**3) should be absent'
    assert series((1/x+1)**3,x,0,3) == 1+x**(-3)+3*x**(-2)+3/x + O(x**3)
    assert series(1/(1+1/x**2),x,0,6) == x**2-x**4+O(x**6, x)

def test_bug2():
    w = Symbol("w")
    p = (w**(-1)+w**(-log(3)*log(2)**(-1)))**(-1)*(3*w**(-log(3)*log(2)**(-1))+2*w**(-1))
    p = p.expand()
    p1 = series(p,w,0,4)
    assert limit(p1,w,0) == 3
    # NO, this does not work
    #assert p1.subs(w, 0) == 3

    # test_exp NO not done by series method
    #p = (1+x)**(1/x)
    #p1 = series(p,0,3)
    #assert p1 == exp(1) - x*exp(1)/2 + O(x**2, x)


def test_exp2():
    pass
    # NO specific to nseries'
    # e = w**(1-log(x)/(log(2) + log(x)))
    # logw = Symbol("logw")
    # e.nseries(w,0,1,logx=logw)

def test_bug3():
    p = (2/x+3/x**2)/(1/x+1/x**2)
    p1 = series(p,x,0,1)
    assert p1 == 3 + O(x)

def test_generalexponent():
    pw=2
    p = (2/x+3/x**pw)/(1/x+1/x**pw)
    p1 = series(p,x,0,1)
    assert p1 == 3 + O(x)
    pw = Rational(1,2)
    p = (2/x+3/x**pw)/(1/x+1/x**pw)
    p1 = series(p,x,0,1)
    assert p1 == 2 + sqrt(x) + O(x)
    p = 1+x**Rational(1,2)
    p1 = series(p,x,0,2)
    assert p1 == 1+x**Rational(1,2)

def test_genexp_x():
    p = 1/(1+x**Rational(1,2))
    p1 = series(p,x,0,2)
    assert p1 == 1+x-x**Rational(1,2)-x**Rational(3,2)+O(x**2, x)

def test_genexp_x2():
    pw = Rational(3,2)
    p = (2/x+3/x**pw)/(1/x+1/x**pw)
    p1 = series(p,x,0,2)
    #p2 = p.series(x,0,2)
    #p2n = 3 - sqrt(x) + x + O(sqrt(x)**3)
    assert p1 == 3 - sqrt(x) + x - x**Rational(3,2) +O(x**2, x)

def test_seriesbug2():
    w = Symbol("w")
    p = ((2*w)/w)**(1+w)
    p1 = series(p,w,0,1)
    assert p1 == 2 + O(w, w)

def test_seriesbug2b():
    w = Symbol("w")
    p = sin(2*w)/w
    w = Symbol("w")
    p1 = series(p,w,0,2)
    assert p1 == 2 + O(w**2, w)

def test_seriesbug2d():
    # use of real=True ?'
    w = Symbol("w", real=True)
    p = log(sin(2*w)/w)
    p1 = series(p,w,0,5)
    assert p1 == log(2) - 2*w**2/3 - 4*w**4/45 + O(w**5)

def test_seriesbug2c():
    w = Symbol("w", real=True)
    p =(sin(2*w)/w)**(1+w)
    assert series(p,w,0,1) == 2 + O(w)
    p1 = series(p,w,0,3)
    p1 = p1.expand()
    assert p1 == 2-Rational(4,3)*w**2+w**2*log(2)**2+2*w*log(2)+O(w**3, w)
    p1 = series(p,w,0,1)
    #p2 = p.series(w,0,1)
    assert p1 == 2 + O(w)
    assert p1.subs(w,0) == 2

def test_expbug4():
    x = Symbol("x", real=True)
    p = (log(sin(2*x)/x)*(1+x))
    p1 = series(p,x,0,2)
    assert p1 == log(2) + x*log(2) + O(x**2, x)

    p = exp(log(sin(2*x)/x)*(1+x))
    p1 = series(p,x,0,2)
    assert p1 == 2 + 2*x*log(2) + O(x**2)

def test_logbug4():
    pass
    # NO series method fails'
    #x = Symbol("x")
    #p = 2 + O(x)
    #p1 = series(p,x,0,2)
    #p2 = p.series(x,0,2)

def test_expbug5():
    p = exp(log(1+x)/x)
    p1 = series(p,x,0,2)
    assert p1 == exp(1) + -exp(1)*x/2 + O(x**2)

def test_sinsinbug():
    p1 = series(sin(sin(x)),x,0,8)
    assert p1 == x-x**3/3+x**5/10-8*x**7/315+O(x**8)

def test_issue159():
    p = x/(exp(x)-1)
    p1 = series(p,x,0,5)
    assert p1 == 1 - x/2 - x**4/720 + x**2/12 + O(x**5)

def test_issue105():
    x = Symbol("x", nonnegative=True)
    p = sin(x**3)**Rational(1,3)
    p1 = series(p,x,0,17)
    assert p1 == x - x**7/18 - x**13/3240 + O(x**17)

    y = Symbol("y")
    p = (1-y**(Rational(1)/2))**(Rational(1)/2)
    p1 = series(p,y,0,2)
    assert p1 == 1 - sqrt(y)/2-y/8-y**Rational(3,2)/16+O(y**2)

    # NO series method fails
    #w = Symbol("w")
    #x = Symbol("x")
    #p = 1/x*(-log(w**(1 + 1/log(3)*log(5))) + log(w + w**(1/log(3)*log(5))))
    #p1 = series(p,w,0,3)

def test_sin():
    x = Symbol("x")
    y = Symbol("y")
    assert series(sin(8*x),x,0,4) == 8*x - 256*x**3/3 + O(x**4)
    assert series(sin(x+y),x,0,1) == sin(y) + O(x)
    assert series(sin(x+y),x,0,2) == sin(y) + cos(y)*x + O(x**2)
    assert series(sin(x+y),x,0,5) == sin(y) + cos(y)*x - sin(y)*x**2/2 - \
                cos(y)*x**3/6 + sin(y)*x**4/24 + O(x**5)

def test_issue416():
    assert series(sin(8*x)/x,x,0,5) == 8 - 256*x**2/3 + 4096*x**4/15 + O(x**5)

def test_issue406():
    p = sin(x)**(-4)*(cos(x)**Rational(1,2)*sin(x)**2 - \
                    cos(x)**Rational(1,3)*sin(x)**2)
    assert series(p,x,0,5) == -Rational(1)/12 - 7*x**2/288 - \
                    43*x**4/10368 + O(x**5)

def test_issue403():
    p = sin(5*x)/sin(2*x)
    assert series(p,x,0,1) == Rational(5,2) + O(x)
    assert series(p,x,0,5) == Rational(5,2) - 35*x**2/4 + 329*x**4/48 + O(x**5)

def test_issue404():
    p = sin(2 + x)/(2 + x)
    p1 = series(p,x,0,2)
    p1 = p1.expand()
    assert p1== sin(2)/2 + x*cos(2)/2 - x*sin(2)/4 + O(x**2)

def test_issue407():
    p = (x + sin(3*x))**(-2)*(x*(x + sin(3*x)) - (x + sin(3*x))*sin(2*x))
    assert series(p,x,0,5) == -Rational(1,4) + 5*x**2/96 + 91*x**4/768 + O(x**5)

def test_issue409():
    x = Symbol("x", real=True)
    p = log(sin(x))
    p1 = series(p,x,0,5)
    assert p1 == log(x) - x**2/6 - x**4/180 + O(x**5)

    p = -log(x) + x*(-log(x) + log(sin(2*x))) + log(sin(2*x))
    p1 = series(p,x,0,5)
    assert p1 == log(2)+log(2)*x-2*x**2/3-2*x**3/3-4*x**4/45+O(x**5)

def test_issue408():
    x = Symbol("x")
    p = x**(-4)*(x**2 - x**2*cos(x)**Rational(1,2))
    p1 = series(p,x,0,5)
    assert p1 == Rational(1,4) + x**2/96 + 19*x**4/5760 + O(x**5)

def test_issue540():
    #taylor collects terms'
    p1 = series(sin(cos(x)),x,0,5)
    assert p1.expand() == sin(1) -x**2*cos(1)/2 - x**4*sin(1)/8 + x**4*cos(1)/24 + O(x**5)

def test_hyperbolic():
    #coth test missing
    x = Symbol("x")
    assert series(sinh(x),x,0,6) == x + x**3/6 + x**5/120 + O(x**6)
    assert series(cosh(x),x,0,5) == 1 + x**2/2 + x**4/24 + O(x**5)
    assert series(tanh(x),x,0,6) == x - x**3/3 + 2*x**5/15 + O(x**6)
    #assert series(cothh(x),x,0,6) == 1/x - x**3/45 + x/3 + 2*x**5/945 + O(x**6)
    assert series(asinh(x),x,0,6) == x - x**3/6 + 3*x**5/40 + O(x**6)
    assert series(acosh(x),x,0,6) == pi*I/2 - I*x - 3*I*x**5/40 - I*x**3/6 + O(x**6)

def test_series2():
    w = Symbol("w", real=True)
    x = Symbol("x", real=True)
    p =  w**(-2)*(w*exp(1/x - w) - w*exp(1/x))
    p1 = series(p,w,0,2)
    assert p1 == -exp(1/x) + w * exp(1/x) / 2  + O(w**2)

def test_series3():
    w = Symbol("w", real=True)
    x = Symbol("x", real=True)
    p = w**(-6)*(w**3*tan(w) - w**3*sin(w))
    p1 = series(p,w,0,2)
    assert p1 == Integer(1)/2 + O(w**2)

def test_bug4():
    w = Symbol("w")
    x = Symbol("x")
    p = x/(w**4 + x**2*w**4 + 2*x*w**4)*w**4
    p1 = series(p,w,0,2)
    #p2 = p.series(w,0,2)
    # the order should be absent
    assert p1 == x/(1 + 2*x + x**2) + O(w**2)

def test_bug5():
    # NO, specific to nseries
    pass
    #w = Symbol("w")
    #x = Symbol("x")
    #l = Symbol('l')
    #e = (-log(w) + log(1 + w*log(x)))**(-2)*w**(-2)*((-log(w) + log(1 + \
    #    x*w))*(-log(w) + log(1 + w*log(x)))*w - x*(-log(w) + log(1 + \
    #        w*log(x)))*w)
    #assert e.nseries(w, n=1, logx=l) == x/w/l + 1/w + O(1, w)
    #assert e.nseries(w, n=2, logx=l) == x/w/l + 1/w - x/l + 1/l*log(x)\
    #        + x*log(x)/l**2 + O(w)

def test_issue1016():
    x = Symbol("x")
    p = sin(x)/(1 - cos(x))
    p1 = series(p,x,0,1)
    assert p1 == 2/x + O(x)
    # O(1/x)  strange notation
    #assert ( sin(x)/(1 - cos(x)) ).nseries(x, n=2) == O(1/x)
    p = sin(x)**2/(1 - cos(x))
    p1 = series(p,x,0,1)
    assert p1 == 2 + O(x)
    # p.series(x,0,1) gives O(x), which is wrong
    #assert ( sin(x)**2/(1 - cos(x)) ).nseries(x, n=2) == O(1, x)

def test_pole():
    x = Symbol("x")
    raises(PoleError, "series(sin(1/x), x, 0, 5)")
    raises(PoleError, "series(sin(1+1/x), x, 0, 5)")
    raises(PoleError, "series(x*sin(1/x), x, 0, 5)")

def test_expsinbug():
    x = Symbol("x")
    assert series(exp(sin(x)),x,0,0) == O(1, x)
    assert series(exp(sin(x)),x,0,1) == 1+O(x)
    assert series(exp(sin(x)),x,0,2) == 1+x+O(x**2)
    assert series(exp(sin(x)),x,0,3) == 1+x+x**2/2+O(x**3)
    assert series(exp(sin(x)),x,0,4) == 1+x+x**2/2+O(x**4)
    assert series(exp(sin(x)),x,0,5) == 1+x+x**2/2-x**4/8+O(x**5)

def test_floor():
    x = Symbol('x')
    assert series(floor(x),x) == 0
    assert series(floor(-x),x) == -1
    assert series(floor(sin(x)),x) == 0
    assert series(floor(sin(-x)),x) == -1
    assert series(floor(x**3),x) == 0
    assert series(floor(-x**3),x) == -1
    assert series(floor(cos(x)),x) == 0
    assert series(floor(cos(-x)),x) == 0
    assert series(floor(5+sin(x)),x) == 5
    assert series(floor(5+sin(-x)),x) == 4
    assert series(floor(x),x,2) == 2
    assert series(floor(-x),x,2) == -3
    x = Symbol('x', negative=True)
    assert series(floor(x+1.5)) == 1

def test_ceiling():
    x = Symbol('x')
    assert series(ceiling(x), x) == 1
    assert series(ceiling(-x), x) == 0
    assert series(ceiling(sin(-x)), x) == 0
    assert series(ceiling(1-cos(x)), x) == 1
    assert series(ceiling(1-cos(-x)), x) == 1
    assert series(ceiling(x), x, 2) == 3
    assert series(ceiling(-x), x, 2) == -2

def test_abs():
    x = Symbol('x')
    a = Symbol('a')
    assert series(abs(x), x, n=4) == x
    assert series(abs(-x), x, n=4) == x
    assert series(abs(x+1), x, n=4) == x+1
    assert series(abs(sin(x)), x, n=4) == x - Rational(1, 6)*x**3 + O(x**4)
    assert series(abs(sin(-x)), x, n=4) == x - Rational(1, 6)*x**3 + O(x**4)
    assert series(abs(x-a), x, 1) ==  Piecewise((x - 1, Eq(1 - a, 0)),
                                            ((x - a)*sign(1 - a), True))

def test_dir():
    x = Symbol('x')
    y = Symbol('y')
    assert series(abs(x), x, 0, dir="+") == x
    assert series(abs(x), x, 0, dir="-") == -x
    assert series(floor(x+2), x, 0, dir="+") == 2
    assert series(floor(x+2), x, 0, dir="-") == 1
    assert series(floor(x+2.2), x, 0, dir="-") == 2
    assert series(ceiling(x+2.2), x, 0, dir="-") == 3
    assert series(sin(x+y), x, 0, dir="-") == series(sin(x+y), x, 0, dir='+')

def test_issue405():
    a = Symbol("a")
    p = asin(a*x)/x
    p1 = series(p,x, 4, n=2)
    assert p1.removeO().subs(x, x - 4) == (
           asin(4*a)/4 -
           (x - 4)*asin(4*a)/16 +
           a*(x - 4)/(4*sqrt(1 - 16*a**2)))

def test_issue1342():
    x, a, b = symbols('x,a,b')
    p = 1/(1+a*x)
    assert series(p, x, 0, 5) == 1 - a*x + a**2*x**2 - a**3*x**3 + \
            a**4*x**4 + O(x**5)
    p = 1/(1+(a+b)*x)
    p1 = series(p, x, 0, 3)
    p1 = p1.expand()
    assert p1 == 1 - a*x - b*x + a**2*x**2 + b**2*x**2 + \
            2*a*b*x**2 + O(x**3)

def test_issue1230():
    assert series(tan(x), x, pi/2, n=3).removeO().subs(x, x - pi/2) == \
           -pi/6 + x/3 - 1/(x - pi/2)
    assert series(cot(x), x, pi, n=3).removeO().subs(x, x - pi) == \
           -x/3 + pi/3 + 1/(x - pi)
    assert limit(tan(x)**tan(2*x), x, pi/4) == exp(-1)

def test_issue2084():
    assert series(abs(x + x**2), n=1) == O(x)
    assert series(abs(x + x**2), n=2) == x + O(x**2)
    # the order should not be there
    assert series((1+x)**2,x, n=6) == 1 + 2*x + x**2 + O(x**6)
    assert series((1 + 1/x)) == 1 + 1/x
    assert Derivative(series(exp(x)), x).doit() == \
           1 + x + x**2/2 + x**3/6 + x**4/24 + O(x**5)
