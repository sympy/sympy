"""Tests for taylor, computing series using lpoly"""
from sympy.polys.ltaylor import taylor, series_reversion, polynomial_degree, poly_truncate, _is_monomial, _get_var_from_term, _factor_var
from sympy.functions.elementary.trigonometric import (cos,sin,atan,acos,acot,tan, asin, cot, atan2)
from sympy.functions.elementary.exponential import (exp,log,LambertW)
from sympy.functions.elementary.hyperbolic import (asinh,acosh,acoth,sinh,cosh,tanh, atanh)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.series.order import O
from sympy.core.numbers import (Rational,E,Integer)
from sympy.core.symbol import Symbol
from sympy import (I,PoleError,floor,ceiling, Piecewise, Eq, sign, symbols, limit, Derivative)
from sympy.core import pi, cacheit, sympify
from sympy.abc import x, y
from sympy.utilities.pytest import raises
from sympy import expand, simplify, Pow, S, Add
from sympy.functions.special.error_functions import erf
from sympy.core.function import Function
from sympy.utilities.pytest import XFAIL

x = Symbol('x')

def test_helper_functions():
    assert _is_monomial(x, x) == (1, 1)
    assert _is_monomial(3*x**2, x) == (2, 3)
    assert _is_monomial(x*(x+1), x) == None
    assert _is_monomial(x**2, x) == (2, 1)
    assert _is_monomial((x+1)**2, x) == None
    assert _is_monomial(2*(x+1), x) == None
    assert _is_monomial(x**x, x) == None
    assert _is_monomial(2*x**x, x) == None
    assert _get_var_from_term((x+1)**2, x) == (0, (x + 1)**2)
    assert _factor_var(1 + x + sin(x)/x**2, x) == (-2, x**2 + x**3 + sin(x))

def test_poly_truncate():
    assert poly_truncate(x.as_poly(x), x, 6) == x
    assert poly_truncate(((1 + x)**4).as_poly(), x, 3) == 6*x**2 + 4*x + 1

def test_polynomial_degree():
    assert polynomial_degree(S.One, x) == 0
    assert polynomial_degree(x, x) == 1
    assert polynomial_degree(-x, x) == 1
    assert polynomial_degree((x+1)*((x**2+1)**2 + 1)**2 + x + 1, x) == 9
    assert polynomial_degree((x+1)*(x**2+1)**2 - x**5, x) == 5
    assert polynomial_degree(1/x, x) == -1
    assert polynomial_degree(x/(1+x), x) == -1
    assert polynomial_degree(sin(x), x) == -1

def test_taylor_QQ1():
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
    assert taylor((S(8)/27 + sin(x))**Rational(1,3), x, 0, 3) == \
      Rational(2,3) + 3*x/4 - 27*x**2/32 + O(x**3)

    assert taylor((1 + 3*x + 2*x**2)**100, x, 0, 5) == \
        1 + 300*x + 44750*x**2 + 4425300*x**3 + 326370825*x**4 + O(x**5)

def test_taylor_QQ2():
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
    assert taylor((1 + x)**Rational(1,4), x, 0, 11) == \
      1 + x/4 - 3*x**2/32 + 7*x**3/128 - 77*x**4/2048 + 231*x**5/8192 - 1463*x**6/65536 + 4807*x**7/262144 - 129789*x**8/8388608 + 447051*x**9/33554432 - 3129357*x**10/268435456 + O(x**11)

    assert taylor((1 + x)**Rational(3,4), x, 0, 11) == \
      1 + 3*x/4 - 3*x**2/32 + 5*x**3/128 - 45*x**4/2048 + 117*x**5/8192 - 663*x**6/65536 + 1989*x**7/262144 - 49725*x**8/8388608 + 160225*x**9/33554432 - 1057485*x**10/268435456 + O(x**11)

    assert taylor((cos(x))**Rational(-3,4), x, 0, 11) == \
      1 + 3*x**2/8 + 17*x**4/128 + 751*x**6/15360 + 63587*x**8/3440640 + 8787601*x**10/1238630400 + O(x**11)

def test_taylor_QQ4():
    assert taylor((exp(x))**Rational(-3,4), x, 0, 7) == \
      1 - 3*x/4 + 9*x**2/32 - 9*x**3/128 + 27*x**4/2048 - 81*x**5/40960 + 81*x**6/327680 + O(x**7)

    assert taylor((1 + log(1+x))**Rational(-3,4), x, 0, 7) == \
      1 - 3*x/4 + 33*x**2/32 - 193*x**3/128 + 4619*x**4/2048 - 139809*x**5/40960 + 5114923*x**6/983040 + O(x**7)

    assert taylor(sqrt((1 + log(1+x))**Rational(-3,4)), x, 0, 7) == \
      1 - 3*x/8 + 57*x**2/128 - 601*x**3/1024 + 26491*x**4/32768 - 1497009*x**5/1310720 + 103245539*x**6/62914560 + O(x**7)

    assert taylor(asinh(x**2), x, 0, 11) == x**2 - x**6/6 + 3*x**10/40 + O(x**11)
    assert taylor(asin(x + x**2), x, 0, 6) == x + x**2 + x**3/6 + x**4/2 + 23*x**5/40 + O(x**6)

def test_taylor_QQ5():
    assert taylor(1/(exp(x)-1), x, 0, 7) == 1/x - S.Half + x/12 - x**3/720 + x**5/30240 + O(x**7)

    assert taylor(sin(x)/sin(tan(x + x**2)), x, 0, 7) == \
      1 - x + 2*x**2/3 - x**3 + 83*x**4/90 - 59*x**5/90 + 1949*x**6/1890 + O(x**7)

def test_taylor_QQ6():
    assert taylor((1 + 1/x)/cos(x), x, 0, 7) == \
      1/x + 1 + x/2 + x**2/2 + 5*x**3/24 + 5*x**4/24 + 61*x**5/720 + 61*x**6/720 + O(x**7)

    assert taylor(cos(x)**sin(x), x, 0, 10) == \
      1 - x**3/2 + x**6/8 - x**7/80 - 37*x**9/1512 + O(x**10)

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
    # series gives O(x**4*log(x)**4) missing terms x**4 * log(x)**k

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

def test_taylor_QQ_root():
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

def test_taylor_QQ_par1():
    assert taylor((x+y)**4, x, 0, 3, pol_pars=[y]) == \
      y**4 + 4*x*y**3 + 6*x**2*y**2 + O(x**3)

    assert taylor((1+x+y)**4, x, 0, 3, pol_pars=[y]) == \
      1 + 4*y + 6*y**2 + 4*y**3 + y**4 + x*(4*y**3 + 12*y**2 + 12*y + 4) + x**2*(6*y**2 + 12*y + 6) + O(x**3)

    p1 = taylor(atan(x*y + x**2), x, 0, 6)
    assert p1 == x*y + x**2 - x**3*y**3/3 - x**4*y**2 + x**5*(y**5/5 - y) + O(x**6)

    assert taylor(sin(x*y), x, 0, 2, pol_pars=[y]) == x*y + O(x**2)
    assert taylor((sin(x) + y)*cos(x), x, 0, 3, pol_pars=[y]) == y + x - x**2*y/2 + O(x**3)

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
    assert p1 == p2

    p1 = taylor(sin(x+40), x, 0, 6)
    assert p1 == sin(40) + x*cos(40) - x**2*sin(40)/2 - x**3*cos(40)/6 + x**4*sin(40)/24 + x**5*cos(40)/120 + O(x**6)

    p1 = taylor(sin(x)/((1 + cos(x)/2)*(2 + cos(x)/2)**S.Half), x, 0, 4)
    assert p1 == 2*sqrt(10)*x/15 + sqrt(10)*x**3/150 + O(x**4)

    assert taylor(1/(cos(x)*(cos(2) + x)), x, 0, 2) == 1/cos(2) - x/cos(2)**2 + O(x**2)

    assert taylor(1/((cos(2) + x)*(1 + cos(x))), x, 0, 2) == 1/(2*cos(2)) - x/(2*cos(2)**2) + O(x**2)

    assert taylor(1/(cos(x/3)*(cos(2) + x/7)), x, 0, 2) == 1/cos(2) - x/(7*cos(2)**2) + O(x**2)

def test_taylor_SR2():
    p1 = taylor(exp(2 - exp(1+x+x**2)), x, 0, 5)
    p2 = exp(-E + 2) - x*exp(-E + 3) + x**2*exp(-E + 4)/2 - 3*x**2*exp(-E + 3)/2 + 3*x**3*exp(-E + 4)/2 - 7*x**3*exp(-E + 3)/6 - x**3*exp(-E + 5)/6 + 55*x**4*exp(-E + 4)/24 + x**4*exp(-E + 6)/24 - 25*x**4*exp(-E + 3)/24 - 3*x**4*exp(-E + 5)/4 + O(x**5)
    # ATT series collects exponentials
    assert p1.expand() == p2.expand()

def test_SR3():
    # example from http://groups.google.com/group/sympy/msg/37fea5ab11302a08
    x,y = symbols('x,y')
    p = (exp(x))/((y/2 + log(2*pi)/2 + x/12 - 1/x - y/x))
    p1 = taylor(taylor(p, x, 0, 3).removeO(), y, 0, 2)
    p2 = -x + x*y - x**2*log(2)/2 - x**2*log(pi)/2 + y*x**2*log(2) + y*x**2*log(pi) - x**2 + y*x**2/2 + O(y**2)

    assert p1 == p2

def test_taylor_SR_logx():
    p1 = taylor(sin(x)**2*log(sin(x) + pi*x)*log(x), x, 0, 6)
    p2 = x**2*log(x)*log(1 + pi) + x**2*log(x)**2 - 120*x**4*log(x)/(720 + 720*pi) - x**4*log(x)*log(1 + pi)/3 - x**4*log(x)**2/3 + O(x**6*log(x)**2)
    assert simplify(p1 - p2) == O(x**6*log(x)**2)

    p1 = taylor(tan(pi*x)**2*log(tan(x)), x, 0, 6)
    p2 = pi**2*x**2*log(x) + x**4*(2*pi**4*log(x)/3 + pi**2/3) + 17*pi**6*x**6*log(x)/45 + O(x**6)
    # series has order O(x**6*log(x)), missing last term
    assert p1 == p2

def test_taylor_SR_root():
    p1 = taylor(sqrt(sin(cos(2)*x)), x, 0, 7)
    assert p1 == sqrt(x)*sqrt(cos(2)) - x**(S(5)/2)*cos(2)**(S(5)/2)/12 + x**(S(9)/2)*cos(2)**(S(9)/2)/1440 - x**(S(13)/2)*cos(2)**(S(13)/2)/24192 + O(x**7)

    p1 = taylor(sqrt(x*sin(cos(2)*x)), x, 0, 7)
    c = cos(2)
    assert p1 == x*sqrt(c) - x**3*c**(S(5)/2)/12 + x**5*c**(S(9)/2)/1440 + O(x**7)

def test_lambert():
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
    p1 = taylor(acos(x), x, 0, 10)
    assert p1 == pi/2 - x - x**3/6 - 3*x**5/40 - 5*x**7/112 - 35*x**9/1152 + O(x**10)

    p1 = taylor(acos(x+x**2), x, 0, 10)
    assert p1 == pi/2 - x - x**2 - x**3/6 - x**4/2 - 23*x**5/40 - 13*x**6/24 - 89*x**7/112 - 17*x**8/16 - 1547*x**9/1152 + O(x**10)

def test_acot():
    p1 = taylor(acot(x+x**2), x, 0, 10)
    p2 = pi/2 - x - x**2 + x**3/3 + x**4 + 4*x**5/5 - 2*x**6/3 - 13*x**7/7 - x**8 + 17*x**9/9 + O(x**10)
    assert p1 == p2

def test_acoth():
    p1 = taylor(acoth(x+x**2), x, 0, 10)
    p2 = I*pi/2 + x + x**2 + x**3/3 + x**4 + 6*x**5/5 + 4*x**6/3 + 15*x**7/7 + 3*x**8 + 37*x**9/9 + O(x**10)
    assert p1 == p2

def test_acosh():
    p1 = taylor(acosh(x), x, 0, 10)
    p2 = I*pi/2 - I*x - I*x**3/6 - 3*I*x**5/40 - 5*I*x**7/112 - 35*I*x**9/1152 + O(x**10)
    assert p1 == p2

def test_sum():
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

@XFAIL
def test_factor_var_from_num_1():
    x = Symbol("x")
    p = ((1 + 1/x)*(1 + 2/x))
    p1 = taylor(((1 + 1/x)*(1 + 2/x)), x, 0, 5)
    # the series method does not give the order, taylor does
    assert p1 == 1 + 3/x + 2/x**2

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
    p1 = taylor(log(sin(2*x)*erf(x)*sqrt(pi)), x, 0, 7, analytic=True)
    assert p1 == 2*log(2) + 2*log(x) - x**2 - 2*x**4/45 - 8*x**6/315 + O(x**7)

    p1 = taylor(sqrt(pi)*erf(sin(x)), x, 0, 10, analytic=True)
    assert p1 == 2*x - x**3 + 11*x**5/20 - 241*x**7/840 + 2777*x**9/20160 + O(x**10)
    p1 = taylor(log(sin(2*x)*erf(x)*sqrt(pi)), x, 0, 10, analytic=True)
    assert p1 == 2*log(2) + 2*log(x) - x**2 - 2*x**4/45 - 8*x**6/315 - 4*x**8/567 + O(x**10)

    p1 = taylor(sqrt(pi)*erf(sin(cos(2)*x)), x, 0, 8, analytic=True)
    assert p1 == 2*x*cos(2) - x**3*cos(2)**3 + 11*x**5*cos(2)**5/20 - 241*x**7*cos(2)**7/840 + O(x**8)

    assert taylor(ftoy1(x), x, 0, 1, analytic=True) == O(x)
    assert taylor(ftoy1(x), x, 0, 2, analytic=True) == x + O(x**2)
    assert taylor(ftoy1(x), x, 0, 4, analytic=True) == x + x**3 + O(x**4)

def test_taylor_series():
    # examples for which currently taylor calls the series method
    p1 = taylor(exp(x*log(x)), x, 0, 3)
    assert p1 == 1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)

    p = taylor(1/cos(x*sqrt(x)), x, 0, 10)
    assert p == 1 + x**3/2 + 5*x**6/24 + 61*x**9/720 + O(x**10)

    p1 = taylor(cos(sqrt(2)*x)*log(tan(x)), x, 0, 6)
    p2 = log(x) + x**2*(-log(x) + S.One/3) + x**4*(log(x)/6 - S(23)/90) - x**6*log(x)/90 + O(x**6)
    assert p1 == p2

    assert taylor(2*exp(1/x), x, 0, 4) == 2*exp(1/x)
    assert taylor(log(1 - cos(x**4)), x, 0, 1) == -log(2) + 8*log(x) + O(x)
    assert taylor(log(x)*log(1 - cos(x**4)), x, 0, 1) == \
            -log(2)*log(x) + 8*log(x)**2 + O(x*log(x))
    assert taylor(log(1 - cos(x**4)), x, 0, 2) == -log(2) + 8*log(x) + O(x**2)
    assert taylor((1 + exp(1/x))**2, x, 0, 2) == exp(2/x) + 2*exp(1/x) + 1
    assert taylor((1 + x)**x * x, x, 0, 2) == x + O(x**2)
    assert taylor((1 + 1/sin(x))/sin(x), x, 0, 2) == x**(-2) + 1/x + S.One/3 + x/6 + O(x**2)
    assert taylor(exp(1/x), x, 0, 2) == exp(1/x)
    assert taylor(exp(1/x)*(1 + sin(x)/x), x, 0, 2) == 2*exp(1/x) + O(x**2*exp(1/x))
    assert taylor(sin((1 + I) + x), x, 0, 2) == sin(1 + I) + x*cos(1 + I) + O(x**2)
    assert taylor(sin(x**(1+x)), x, 0, 2) == exp((x + 1)*log(x)) + O(x**2)
    assert taylor(atan2(x, y), x, 0, 2, analytic=True) == atan2(0, y) + x/y + O(x**2)
    assert taylor(erf(x + 1), x, 0, 2, analytic=True) == erf(1) + 2*x*exp(-1)/sqrt(pi) + O(x**2)
    assert taylor(cot(x), x, 0, 4, analytic=True) == 1/x - x/3 - x**3/45 + O(x**4)
    assert taylor(ftoy2(x), x, 0, 2, analytic=True) == 2/x**3 + 2/x**2 + 2/x + O(x**2)
    assert taylor(ftoy3(x), x, 0, 2, analytic=True) == 1 + x + O(x**2)

@XFAIL
def test_taylor_series_fail1():
    # taylor calls the series method, which fails
    assert taylor(log(1 - cos(x**4)), x, 0, 3) == -log(2) + 8*log(x) + O(x**3)

def test_rdeco():
    assert taylor(log(1 - cos(x**4)), x, 0, 3, rdeco=3) == -log(2) + 8*log(x) + O(x**3)

#######################################################
# tests from test_nseries.py, adapted

def test_nseries():
    """tests adapted from test_nseries.py"""
    assert taylor(x, x, n=5) == x
    assert taylor(y, x, n=5) == y

def test_simple_1():
    assert taylor(x, x, n=5) == x
    y = Symbol("y")
    assert taylor(y, x, n=5) == y
    assert taylor(1/(x*y),x, n=5) == 1/(x*y)
    assert taylor(Rational(3, 4), x, n=5) == Rational(3,4)
    assert taylor(x) == x

def test_mul_0():
    assert taylor(x*log(x), x, n=5) == x*log(x)

def test_mul_1():
    assert taylor(x*log(2+x), x, n=4) == x*log(2) + x**2/2 - x**3/8 + O(x**4)
    assert taylor(x*log(1+x), x, n=4) == x**2 - x**3/2 + O(x**4)

def test_pow_0():
    assert taylor(x**2, x, n=5) == x**2
    assert taylor(1/x, x, n=5) == 1/x
    assert taylor(1/x**2, x, n=5) == 1/x**2
    assert taylor(x**(Rational(2,3)), x, n=5) == (x**(Rational(2,3)))
    assert taylor(sqrt(x)**3, x, n=5) == sqrt(x)**3

def test_pow_1():
    assert taylor((1+x)**2, x, n=5) == 1+2*x+x**2

def test_geometric_1():
    assert taylor(1/(1-x), x, n=5) == 1+x+x**2+x**3+x**4+O(x**5)
    assert taylor(x/(1-x), x, n=5) == x + x**2 + x**3 + x**4 + O(x**5)
    assert taylor(x**3/(1-x), x, n=5) == x**3 + x**4 + O(x**5)

def test_sqrt_1():
    assert taylor(sqrt(1+x), x, n=5) == 1+x/2-x**2/8+x**3/16-5*x**4/128+O(x**5)

def test_exp_1():
    assert taylor(exp(x), x, n=5) == 1+x+x**2/2+x**3/6+x**4/24 + O(x**5)
    assert taylor(exp(x), x, n=12) == 1+x+x**2/2+x**3/6+x**4/24+x**5/120+  \
               x**6/720+x**7/5040+x**8/40320+x**9/362880+x**10/3628800+  \
               x**11/39916800 + O(x**12)
    assert taylor(exp(1/x), x, n=5) == exp(1/x)
    assert taylor(exp(1/(1+x)), x, n=4) == \
            E - E*x + 3*E*x**2/2 - 13*E*x**3/6 + O(x**4)
    assert taylor(exp(2+x), x, n=5) == \
      exp(2) + x*exp(2) + x**2*exp(2)/2 + x**3*exp(2)/6 + x**4*exp(2)/24 + O(x**5)

def test_exp_sqrt_1():
    p1 = taylor(exp(1+sqrt(x)), x, n=3)
    s = sqrt(x)
    p2 = E + E*x/2 + E*s + E*x**2/24 + E*s**3/6 + E*s**5/120 + O(x**3)
    assert p1 == p2

def test_power_x_x1():
    assert taylor(exp(x*log(x)), x, n=4) == \
      1+x*log(x)+x**2*log(x)**2/2+x**3*log(x)**3/6 + O(x**4*log(x)**4)

def test_power_x_x2():
    assert taylor((x**x), x, n=4) == \
       1 + x*log(x) + x**2*log(x)**2/2 + x**3*log(x)**3/6 + O(x**4*log(x)**4)

def test_log_singular1():
    assert taylor(log(1+1/x), x, n=5) == \
      x - log(x) - x**2/2 + x**3/3 - x**4/4 + O(x**5)

def test_log_power1():
    assert taylor(1 / (1/x + x ** (log(3)/log(2))), x, n=5) == \
      x - x**(2 + log(3)/log(2)) + O(x**5)

#def test_log_series():
#    l = Symbol('l')
#    e = 1/(1-log(x))
#    assert e.nseries(x, n=5, logx=l) == 1/(1-l)
#    (NO) not done since the series method does not do it

def test_log2():
    assert taylor(log(-1/x), x, n=5) == -log(x) + log(-1)

#def test_log3():
#    l = Symbol('l')
#    e = 1/log(-1/x)
#    assert e.nseries(x, n=4, logx=l) == 1/(-l + log(-1))
#    (NO) not done since the series method does not do it

def test_series1():
    x = Symbol("x")
    p = sin(x)
    assert taylor(p, x, 0, 0) != 0
    assert taylor(p, x, 0, 0) == O(1,x)
    assert taylor(p, x, 0, 1) == O(x,x)
    assert taylor(p, x, 0, 2) == x + O(x**2, x)
    assert taylor(p, x, 0, 3) == x + O(x**3, x)
    assert taylor(p, x, 0, 4) == x-x**3/6 + O(x**4, x)
    p = (exp(x)-1)/x
    #p2n = 1+x/2+O(x**2, x)
    assert taylor(p, x, 0, 3) == 1 + x/2 + x**2/6 + O(x**3, x)
    assert taylor(x, x, 0, 3) == x
    # in test_series1_failing
    assert taylor(x, x, 0, 0) == O(1)  # (NO) O(1, x)
    assert taylor(x, x, 0, 1) == O(x, x)

def test_seriesbug1():
    assert taylor(1/x, x, 0, 3) == 1/x
    assert taylor(x + 1/x, x, 0, 3) == x + 1/x

def test_series2x():
    assert taylor((x+1)**(-2), x, 0, 4) == 1-2*x+3*x**2-4*x**3+O(x**4, x)
    assert taylor((x+1)**(-1), x, 0, 4) == 1-x+x**2-x**3+O(x**4, x)
    assert taylor((x+1)**0, x, 0, 3) == 1
    assert taylor((x+1)**1, x, 0, 3) == 1+x
    assert taylor((x+1)**2, x, 0, 3) == 1+2*x+x**2
    #p2n = 1 + 3*x + 3*x**2 + x**3; with series metod there in O(x**3)
    assert taylor((x+1)**3, x, 0, 3) == 1+3*x+3*x**2 + O(x**3)
    assert taylor(1/(1+x), x, 0, 4) == 1-x+x**2-x**3+O(x**4, x)
    assert taylor(x+3/(1+2*x), x, 0, 4) == 3-5*x+12*x**2-24*x**3+O(x**4, x)
    # with series metod there in O(x**3)
    assert taylor((1/x+1)**3, x, 0, 3) == 1+x**(-3)+3*x**(-2)+3/x + O(x**3)
    assert taylor(1/(1+1/x**2),x, 0, 6) == x**2-x**4+O(x**6, x)

def test_bug2():
    w = Symbol("w")
    p = (w**(-1)+w**(-log(3)*log(2)**(-1)))**(-1)*(3*w**(-log(3)*log(2)**(-1))+2*w**(-1))
    p = p.expand()
    p1 = taylor(p,w,0,4)
    assert limit(p1,w,0) == 3
    # NO, this does not work
    #assert p1.subs(w, 0) == 3

def test_exp():
    p = (1+x)**(1/x)
    assert taylor(p, x, n=3) == E - E*x/2 + 11*E*x**2/24 + O(x**3)

#def test_exp2():
# NO specific to nseries
    # w = Symbol("w")
    # e = w**(1-log(x)/(log(2) + log(x)))
    # logw = Symbol("logw")
    # e.nseries(w,0,1,logx=logw)

def test_bug3():
    p = (2/x+3/x**2)/(1/x+1/x**2)
    assert taylor(p,x,0,1) == 3 + O(x)

def test_generalexponent():
    pw=2
    p = (2/x+3/x**pw)/(1/x+1/x**pw)
    assert taylor(p, x, 0, 1) == 3 + O(x)
    pw = Rational(1,2)
    assert taylor((2/x+3/x**pw)/(1/x+1/x**pw), x, 0, 1) == 2 + sqrt(x) + O(x)
    assert taylor(1+x**pw, x, 0, 2) == 1+x**pw
    assert taylor(1+sqrt(x), x, 0, 4) == 1+sqrt(x)

def test_genexp_x():
    p = 1/(1+x**Rational(1,2))
    p1 = taylor(p, x, 0, 2)
    assert p1 == 1+x-x**Rational(1,2)-x**Rational(3,2)+O(x**2, x)

def test_genexp_x2():
    pw = Rational(3,2)
    p = (2/x+3/x**pw)/(1/x+1/x**pw)
    p1 = taylor(p, x, 0, 2)
    assert p1 == 3 - sqrt(x) + x - x**Rational(3,2) + O(x**2, x)

def test_seriesbug2():
    w = Symbol("w")
    p = ((2*w)/w)**(1+w)
    assert taylor(p, w, 0, 1) == 2 + O(w, w)

def test_seriesbug2b():
    w = Symbol("w")
    p = sin(2*w)/w
    assert taylor(p, w, 0, 2) == 2 + O(w**2, w)

def test_seriesbug2d():
    # use of real=True ?
    w = Symbol("w", real=True)
    p = log(sin(2*w)/w)
    assert taylor(p, w, n=5) == log(2) - 2*w**2/3 - 4*w**4/45 + O(w**5)

def test_seriesbug2c():
    w = Symbol("w", real=True)
    p = (sin(2*w)/w)**(1+w)
    assert taylor(p, w, 0, 1) == 2 + O(w)
    p1 = taylor(p, w, 0, 3)
    assert p1 == 2 + 2*w*log(2) + w**2*(-S(4)/3 + log(2)**2) + O(w**3)

    p1 = taylor(p, w, 0, 1)
    assert p1 == 2 + O(w)
    assert p1.subs(w, 0) == 2

def test_expbug4():
    x = Symbol("x", real=True)
    p = (log(sin(2*x)/x)*(1+x))
    p1 = taylor(p,x,0,2)
    assert p1 == log(2) + x*log(2) + O(x**2, x)

    p = exp(log(sin(2*x)/x)*(1+x))
    p1 = taylor(p,x,0,2)
    assert p1 == 2 + 2*x*log(2) + O(x**2)

@XFAIL
def test_expbug4_failing_1():
    x = Symbol("x", real=True)
    assert taylor(exp(log(2)+O(x)), x, 0, 2) == 2 + O(x**2, x)

@XFAIL
def test_expbug4_failing_2():
    x = Symbol("x", real=True)
    assert taylor(((2+O(x))**(1+x)), x, 0, 2) == 2 + O(x**2, x)

@XFAIL
def test_expbug5_failing():
    x = Symbol("x", real=True)
    assert taylor(exp(O(x)), x, 0, 2) == 1 + O(x**2, x)

def test_logbug4():
    assert taylor(log(sin(2*x)/x)*(1+x), x, 0, 2) == log(2) + x*log(2) + O(x**2)
    assert taylor(exp(log(sin(2*x)/x)*(1+x)), x, 0, 2) == 2 + 2*x*log(2) + O(x**2)

def test_expbug5():
    assert taylor(exp(log(1+x)/x), x, 0, 2) == exp(1) + -exp(1)*x/2 + O(x**2)

def test_sinsinbug():
    assert taylor(sin(sin(x)), x, 0, 8) == x-x**3/3+x**5/10-8*x**7/315+O(x**8)

def test_issue159():
    assert taylor(x/(exp(x)-1), x, 0, 5) == 1 - x/2 - x**4/720 + x**2/12 + O(x**5)

def test_issue105():
    x = Symbol("x", nonnegative=True)
    p = sin(x**3)**Rational(1,3)
    assert taylor(p,x,0,17) == x - x**7/18 - x**13/3240 + O(x**17)

def test_issue125():
    y = Symbol("y")
    assert taylor(sqrt(1-sqrt(y)), y, 0, 2) == 1 - sqrt(y)/2-y/8-sqrt(y)**3/16+O(y**2)

def test_issue364():
    w,x,i = symbols('w,x,i')
    r = log(5)/log(3)
    p = w**(-1 + r)
    e = 1/x*(-log(w**(1 + r)) + log(w + w**r))
    e_ser = -r*log(w)/x + p/x - p**2/(2*x) + O(p**3)
    assert e.nseries(w, n=3) == e_ser

def test_sin():
    y = Symbol("y")
    assert taylor(sin(8*x), x, n=4) == 8*x - 256*x**3/3 + O(x**4)
    assert taylor(sin(x+y), x, n=1) == sin(y) + O(x)
    assert taylor(sin(x+y), x, n=2) == sin(y) + cos(y)*x + O(x**2)
    assert taylor(sin(x+y), x, n=5) == sin(y) + cos(y)*x - sin(y)*x**2/2 - \
                cos(y)*x**3/6 + sin(y)*x**4/24 + O(x**5)

def test_issue416():
    assert taylor(sin(8*x)/x, x, n=5) == 8 - 256*x**2/3 + 4096*x**4/15 + O(x**5)

def test_issue406():
    p = sin(x)**(-4)*(sqrt(cos(x))*sin(x)**2 - \
            cos(x)**Rational(1,3)*sin(x)**2)
    assert taylor(p,x,0,5) == -Rational(1)/12 - 7*x**2/288 - \
                    43*x**4/10368 + O(x**5)

def test_issue402():
    a = Symbol("a")
    e = x**(-2)*(x*sin(a + x) - x*sin(a))
    assert taylor(e, x, n=4) == cos(a) - sin(a)*x/2 - cos(a)*x**2/6 + \
            sin(a)*x**3/24 + O(x**4)
    e = x**(-2)*(x*cos(a + x) - x*cos(a))
    assert taylor(e, x, n=4) == -sin(a) - cos(a)*x/2 + sin(a)*x**2/6 + \
            cos(a)*x**3/24 + O(x**4)

def test_issue403():
    p = sin(5*x)/sin(2*x)
    assert taylor(p,x, n=1) == Rational(5,2) + O(x)
    assert taylor(p, x, n=5) == Rational(5,2) - 35*x**2/4 + 329*x**4/48 + O(x**5)

def test_issue404():
    p = sin(2 + x)/(2 + x)
    p1 = taylor(p, x, n=2)
    assert p1 == sin(2)/2 + x*(-sin(2)/4 + cos(2)/2) + O(x**2)

def test_issue407():
    p = (x + sin(3*x))**(-2)*(x*(x + sin(3*x)) - (x + sin(3*x))*sin(2*x))
    assert taylor(p, x, n=5) == -Rational(1,4) + 5*x**2/96 + 91*x**4/768 + O(x**5)

def test_issue409():
    x = Symbol("x", real=True)
    assert taylor(log(sin(x)), x, n=5) == log(x) - x**2/6 - x**4/180 + O(x**5)
    p = -log(x) + x*(-log(x) + log(sin(2*x))) + log(sin(2*x))
    assert taylor(p, x, n=5) == log(2)+log(2)*x-2*x**2/3-2*x**3/3-4*x**4/45+O(x**5)

def test_issue408():
    p = x**(-4)*(x**2 - x**2*sqrt(cos(x)))
    assert taylor(p, x, n=5) == Rational(1,4) + x**2/96 + 19*x**4/5760 + O(x**5)

def test_issue540():
    #taylor collects terms'
    p1 = taylor(sin(cos(x)), x, n=5)
    assert p1 == sin(1) - x**2*cos(1)/2 + x**4*(-sin(1)/8 + cos(1)/24) + O(x**5)

def test_hyperbolic():
    #coth test missing
    assert taylor(sinh(x), x, n=6) == x + x**3/6 + x**5/120 + O(x**6)
    assert taylor(cosh(x), x, n=5) == 1 + x**2/2 + x**4/24 + O(x**5)
    assert taylor(tanh(x), x, n=6) == x - x**3/3 + 2*x**5/15 + O(x**6)
    #assert series(cothh(x), x, n=6) == 1/x - x**3/45 + x/3 + 2*x**5/945 + O(x**6)
    assert taylor(asinh(x), x, n=6) == x - x**3/6 + 3*x**5/40 + O(x**6)
    assert taylor(acosh(x), x, n=6) == pi*I/2 - I*x - 3*I*x**5/40 - I*x**3/6 + O(x**6)

def test_series2():
    w = Symbol("w", real=True)
    x = Symbol("x", real=True)
    p =  w**(-2)*(w*exp(1/x - w) - w*exp(1/x))
    assert taylor(p, w, n=2) == -exp(1/x) + w * exp(1/x)/2  + O(w**2)

def test_series3():
    w = Symbol("w", real=True)
    p = w**(-6)*(w**3*tan(w) - w**3*sin(w))
    assert taylor(p, w, n=2) == Integer(1)/2 + O(w**2)

def test_bug4():
    w = Symbol("w")
    p = x/(w**4 + x**2*w**4 + 2*x*w**4)*w**4
    # the order should be absent
    assert taylor(p, w, n=2) == x/(1 + 2*x + x**2) + O(w**2)

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
    assert taylor(sin(x)/(1 - cos(x)), x, n=1) == 2/x + O(x)
    assert taylor(sin(x)**2/(1 - cos(x)), x, n=1) == 2 + O(x)

def test_pole():
    raises(PoleError, "taylor(sin(1/x), x, 0, 5)")
    raises(PoleError, "taylor(sin(1+1/x), x, 0, 5)")
    raises(PoleError, "taylor(x*sin(1/x), x, 0, 5)")

def test_expsinbug():
    assert taylor(exp(sin(x)), x, 0, 0) == O(1, x)
    assert taylor(exp(sin(x)), x, 0, 1) == 1+O(x)
    assert taylor(exp(sin(x)), x, 0, 2) == 1+x+O(x**2)
    assert taylor(exp(sin(x)), x, 0, 3) == 1+x+x**2/2+O(x**3)
    assert taylor(exp(sin(x)), x, 0, 4) == 1+x+x**2/2+O(x**4)
    assert taylor(exp(sin(x)), x, 0, 5) == 1+x+x**2/2-x**4/8+O(x**5)

def test_floor():
    x = Symbol('x')
    assert taylor(floor(x), x) == 0
    assert taylor(floor(-x), x) == -1
    assert taylor(floor(sin(x)), x) == 0
    assert taylor(floor(sin(-x)), x) == -1
    assert taylor(floor(x**3), x) == 0
    assert taylor(floor(-x**3), x) == -1
    assert taylor(floor(cos(x)), x) == 0
    assert taylor(floor(cos(-x)), x) == 0
    assert taylor(floor(5+sin(x)), x) == 5
    assert taylor(floor(5+sin(-x)), x) == 4
    assert taylor(floor(x), x ,2) == 2
    assert taylor(floor(-x), x, 2) == -3

    x = Symbol('x', negative=True)
    assert taylor(floor(x+1.5)) == 1

def test_ceiling():
    x = Symbol('x')
    assert taylor(ceiling(x), x) == 1
    assert taylor(ceiling(-x), x) == 0
    assert taylor(ceiling(sin(-x)), x) == 0
    assert taylor(ceiling(1-cos(x)), x) == 1
    assert taylor(ceiling(1-cos(-x)), x) == 1
    assert taylor(ceiling(x), x, 2) == 3
    assert taylor(ceiling(-x), x, 2) == -2

def test_abs():
    x = Symbol('x')
    a = Symbol('a')
    assert taylor(abs(x), x, n=4) == x
    assert taylor(abs(-x), x, n=4) == x
    assert taylor(abs(x+1), x, n=4) == x+1
    assert taylor(abs(sin(x)), n=4) == x - Rational(1, 6)*x**3 + O(x**4)
    assert taylor(abs(sin(-x)), n=4) == x - Rational(1, 6)*x**3 + O(x**4)
    assert taylor(abs(x-a), x, 1) ==  Piecewise((x - 1, Eq(1 - a, 0)),
                                            ((x - a)*sign(1 - a), True))
def test_dir():
    x = Symbol('x')
    y = Symbol('y')
    assert taylor(abs(x), x, 0, dir="+") == x
    assert taylor(abs(x), x, 0, dir="-") == -x
    assert taylor(floor(x+2), x, 0, dir="+") == 2
    assert taylor(floor(x+2), x, 0, dir="-") == 1
    assert taylor(floor(x+2.2), x, 0, dir="-") == 2
    assert taylor(ceiling(x+2.2), x, 0, dir="-") == 3
    assert taylor(sin(x+y), x, 0, dir="-") == taylor(sin(x+y), x, 0, dir='+')

def test_issue405():
    a = Symbol("a")
    p = asin(a*x)/x
    p1 = taylor(p,x, 4, 2)
    assert p1.removeO().subs(x, x - 4) == (
           asin(4*a)/4 -
           (x - 4)*asin(4*a)/16 +
           a*(x - 4)/(4*sqrt(1 - 16*a**2)))

def test_issue1342():
    x, a, b = symbols('x,a,b')
    p = 1/(1+a*x)
    assert taylor(p, x, 0, 5) == 1 - a*x + a**2*x**2 - a**3*x**3 + \
            a**4*x**4 + O(x**5)
    p = 1/(1+(a+b)*x)
    p1 = taylor(p, x, 0, 3)
    p2 = 1 + x*(-a - b) + x**2*(a**2 + 2*a*b + b**2) + O(x**3)
    assert p1 == p2

def test_issue1230():
    assert taylor(tan(x), x, pi/2, 3).removeO().subs(x, x - pi/2) == \
           -pi/6 + x/3 - 1/(x - pi/2)
    assert taylor(cot(x), x, pi, 3).removeO().subs(x, x - pi) == \
           -x/3 + pi/3 + 1/(x - pi)
    assert limit(tan(x)**tan(2*x), x, pi/4) == exp(-1)

def test_issue2084():
    assert taylor(abs(x + x**2), n=1) == O(x)
    assert taylor(abs(x + x**2), n=2) == x + O(x**2)
    assert taylor((1+x)**2, x, n=6) == 1 + 2*x + x**2
    assert taylor((1 + 1/x)) == 1 + 1/x
    assert Derivative(taylor(exp(x)), x).doit() == 1 + x + x**2/2 + x**3/6 + x**4/24 + O(x**5)

def test_series_reversion():
    from sympy.abc import x, y, z
    from sympy.series.order import O
    p = taylor(sin(x), x, 0, 4)
    assert series_reversion(p, [x, y]) == taylor(asin(y), y, 0, 4)
    p = x + x**2/2 + x**2*y + x**3/3 + x**3*y + O(x**4)
    p1 = series_reversion(p, [x, y, z])
    assert p.removeO().subs(x,p1).series(z, 0, 4) == z + O(z**4)
