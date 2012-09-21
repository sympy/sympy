from sympy import (symbols, Symbol, nan, oo, zoo, I, sinh, sin, acot, pi, atan,
        acos, Rational, sqrt, asin, acot, cot, coth, E, S, tan, tanh, cos,
        cosh, atan2, exp, log, asinh, acoth, atanh, O, cancel, Matrix, re, im,
        Float,Pow)

from sympy.utilities.pytest import XFAIL

def test_sin():
    x, y = symbols('x,y')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert sin(nan) == nan

    assert sin(oo*I) == oo*I
    assert sin(-oo*I) == -oo*I
    assert sin(oo).args[0] == oo

    assert sin(0) == 0

    assert sin(asin(x)) == x
    assert sin(atan(x)) == x / sqrt(1 + x**2)
    assert sin(acos(x)) == sqrt(1 - x**2)
    assert sin(acot(x)) == 1 / (sqrt(1 + 1 / x**2) * x)
    assert sin(atan2(y, x)) == y / sqrt(x**2 + y**2)

    assert sin(pi*I) == sinh(pi)*I
    assert sin(-pi*I) == -sinh(pi)*I
    assert sin(-2*I) == -sinh(2)*I

    assert sin(pi) == 0
    assert sin(-pi) == 0
    assert sin(2*pi) == 0
    assert sin(-2*pi) == 0
    assert sin(-3*10**73*pi) == 0
    assert sin(7*10**103*pi) == 0

    assert sin(pi/2) == 1
    assert sin(-pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(7*pi/2) == -1

    assert sin(pi/3) == S.Half*sqrt(3)
    assert sin(-2*pi/3) == -S.Half*sqrt(3)

    assert sin(pi/4) == S.Half*sqrt(2)
    assert sin(-pi/4) == -S.Half*sqrt(2)
    assert sin(17*pi/4) == S.Half*sqrt(2)
    assert sin(-3*pi/4) == -S.Half*sqrt(2)

    assert sin(pi/6) == S.Half
    assert sin(-pi/6) == -S.Half
    assert sin(7*pi/6) == -S.Half
    assert sin(-5*pi/6) == -S.Half

    assert sin(1*pi/5) == sqrt((5 - sqrt(5)) / 8)
    assert sin(2*pi/5) == sqrt((5 + sqrt(5)) / 8)
    assert sin(3*pi/5) == sin(2*pi/5)
    assert sin(4*pi/5) == sin(1*pi/5)
    assert sin(6*pi/5) == -sin(1*pi/5)
    assert sin(8*pi/5) == -sin(2*pi/5)

    assert sin(-1273*pi/5) == -sin(2*pi/5)

    assert sin(104*pi/105) == sin(pi/105)
    assert sin(106*pi/105) == -sin(pi/105)

    assert sin(-104*pi/105) == -sin(pi/105)
    assert sin(-106*pi/105) == sin(pi/105)

    assert sin(x*I) == sinh(x)*I

    assert sin(k*pi) == 0
    assert sin(17*k*pi) == 0

    assert sin(k*pi*I) == sinh(k*pi)*I

    assert sin(r).is_real == True

    assert isinstance(sin( re(x) - im(y)), sin) == True
    assert isinstance(sin(-re(x) + im(y)), sin) == False

def test_sin_series():
    x = Symbol('x')
    assert sin(x).series(x, 0, 9) == \
                    x - x**3/6 + x**5/120 - x**7/5040 + O(x**9)

def test_sin_rewrite():
    x = Symbol('x')
    assert sin(x).rewrite(exp) == -I*(exp(I*x) - exp(-I*x))/2
    assert sin(x).rewrite(tan) == 2*tan(x/2)/(1 + tan(x/2)**2)
    assert sin(x).rewrite(cot) == 2*cot(x/2)/(1 + cot(x/2)**2)
    assert sin(sinh(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, sinh(3)).n()
    assert sin(cosh(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, cosh(3)).n()
    assert sin(tanh(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, tanh(3)).n()
    assert sin(coth(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, coth(3)).n()
    assert sin(sin(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, sin(3)).n()
    assert sin(cos(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, cos(3)).n()
    assert sin(tan(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, tan(3)).n()
    assert sin(cot(x)).rewrite(exp).subs(x, 3).n() == sin(x).rewrite(exp).subs(x, cot(3)).n()
    assert sin(log(x)).rewrite(Pow)  == I*x**-I / 2 - I*x**I /2

def test_sin_expansion():
    # Note: these formulas are not unique.  The ones here come from the
    # Chebyshev formulas.
    x, y = symbols('x y')
    assert sin(x + y).expand(trig=True) == sin(x)*cos(y) + cos(x)*sin(y)
    assert sin(x - y).expand(trig=True) == sin(x)*cos(y) - cos(x)*sin(y)
    assert sin(y - x).expand(trig=True) == cos(x)*sin(y) - sin(x)*cos(y)
    assert sin(2*x).expand(trig=True) == 2*sin(x)*cos(x)
    assert sin(3*x).expand(trig=True) == -4*sin(x)**3 + 3*sin(x)
    assert sin(4*x).expand(trig=True) == -8*sin(x)**3*cos(x) + 4*sin(x)*cos(x)
    assert sin(2).expand(trig=True) == 2*sin(1)*cos(1)
    assert sin(3).expand(trig=True) == -4*sin(1)**3 + 3*sin(1)

def test_trig_symmetry():
    x = Symbol('x')
    y = Symbol('y')
    k = Symbol('k', integer=True)

    assert sin(-x) == -sin(x)
    assert cos(-x) == cos(x)
    assert tan(-x) == -tan(x)
    assert cot(-x) == -cot(x)
    assert sin(x+pi) == -sin(x)
    assert sin(x+2*pi) == sin(x)
    assert sin(x+3*pi) == -sin(x)
    assert sin(x+4*pi) == sin(x)
    assert sin(x-5*pi) == -sin(x)
    assert cos(x+pi) == -cos(x)
    assert cos(x+2*pi) == cos(x)
    assert cos(x+3*pi) == -cos(x)
    assert cos(x+4*pi) == cos(x)
    assert cos(x-5*pi) == -cos(x)
    assert tan(x+pi) == tan(x)
    assert tan(x-3*pi) == tan(x)
    assert cot(x+pi) == cot(x)
    assert cot(x-3*pi) == cot(x)
    assert sin(pi/2-x) == cos(x)
    assert sin(3*pi/2-x) == -cos(x)
    assert sin(5*pi/2-x) == cos(x)
    assert cos(pi/2-x) == sin(x)
    assert cos(3*pi/2-x) == -sin(x)
    assert cos(5*pi/2-x) == sin(x)
    assert tan(pi/2-x) == cot(x)
    assert tan(3*pi/2-x) == cot(x)
    assert tan(5*pi/2-x) == cot(x)
    assert cot(pi/2-x) == tan(x)
    assert cot(3*pi/2-x) == tan(x)
    assert cot(5*pi/2-x) == tan(x)
    assert sin(pi/2+x) == cos(x)
    assert cos(pi/2+x) == -sin(x)
    assert tan(pi/2+x) == -cot(x)
    assert cot(pi/2+x) == -tan(x)

def test_cos():
    x, y = symbols('x,y')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert cos(nan) == nan

    assert cos(oo*I) == oo
    assert cos(-oo*I) == oo

    assert cos(0) == 1

    assert cos(acos(x)) == x
    assert cos(atan(x)) == 1 / sqrt(1 + x**2)
    assert cos(asin(x)) == sqrt(1 - x**2)
    assert cos(acot(x)) == 1 / sqrt(1 + 1 / x**2)
    assert cos(atan2(y, x)) == x / sqrt(x**2 + y**2)

    assert cos(pi*I) == cosh(pi)
    assert cos(-pi*I) == cosh(pi)
    assert cos(-2*I) == cosh(2)

    assert cos(pi/2) == 0
    assert cos(-pi/2) == 0
    assert cos(pi/2) == 0
    assert cos(-pi/2) == 0
    assert cos((-3*10**73+1)*pi/2) == 0
    assert cos((7*10**103+1)*pi/2) == 0

    assert cos(pi) == -1
    assert cos(-pi) == -1
    assert cos(2*pi)==1
    assert cos(5*pi) == -1
    assert cos(8*pi) == 1

    assert cos(pi/3) == S.Half
    assert cos(-2*pi/3) == -S.Half

    assert cos(pi/4) == S.Half*sqrt(2)
    assert cos(-pi/4) == S.Half*sqrt(2)
    assert cos(11*pi/4) == -S.Half*sqrt(2)
    assert cos(-3*pi/4) == -S.Half*sqrt(2)

    assert cos(pi/6) == S.Half*sqrt(3)
    assert cos(-pi/6) == S.Half*sqrt(3)
    assert cos(7*pi/6) == -S.Half*sqrt(3)
    assert cos(-5*pi/6) == -S.Half*sqrt(3)

    assert cos(1*pi/5) == (sqrt(5) + 1)/4
    assert cos(2*pi/5) == (sqrt(5) - 1)/4
    assert cos(3*pi/5) == -cos(2*pi/5)
    assert cos(4*pi/5) == -cos(1*pi/5)
    assert cos(6*pi/5) == -cos(1*pi/5)
    assert cos(8*pi/5) == cos(2*pi/5)

    assert cos(-1273*pi/5) == -cos(2*pi/5)

    assert cos(104*pi/105) == -cos(pi/105)
    assert cos(106*pi/105) == -cos(pi/105)

    assert cos(-104*pi/105) == -cos(pi/105)
    assert cos(-106*pi/105) == -cos(pi/105)

    assert cos(x*I) == cosh(x)
    assert cos(k*pi*I) == cosh(k*pi)

    assert cos(r).is_real == True

    assert cos(k*pi) == (-1)**k
    assert cos(2*k*pi) == 1

def test_issue_3091():
    c = Float('123456789012345678901234567890.25', '')
    for cls in [sin, cos, tan, cot]:
        assert cls(c*pi) == cls(pi/4)
        assert cls(4.125*pi) == cls(pi/8)
        assert cls(4.7*pi) == cls((4.7 % 2)*pi)

def test_cos_series():
    x = Symbol('x')
    assert cos(x).series(x, 0, 9) == \
                    1 - x**2/2 + x**4/24 - x**6/720 + x**8/40320 + O(x**9)

def test_cos_rewrite():
    x = Symbol('x')
    assert cos(x).rewrite(exp) == exp(I*x)/2 + exp(-I*x)/2
    assert cos(x).rewrite(tan) == (1 - tan(x/2)**2)/(1 + tan(x/2)**2)
    assert cos(x).rewrite(cot) == -(1 - cot(x/2)**2)/(1 + cot(x/2)**2)
    assert cos(sinh(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, sinh(3)).n()
    assert cos(cosh(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, cosh(3)).n()
    assert cos(tanh(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, tanh(3)).n()
    assert cos(coth(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, coth(3)).n()
    assert cos(sin(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, sin(3)).n()
    assert cos(cos(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, cos(3)).n()
    assert cos(tan(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, tan(3)).n()
    assert cos(cot(x)).rewrite(exp).subs(x, 3).n() == cos(x).rewrite(exp).subs(x, cot(3)).n()
    assert cos(log(x)).rewrite(Pow) == x**I/2 + x**-I/2

def test_cos_expansion():
    x, y = symbols('x y')
    assert cos(x + y).expand(trig=True) == cos(x)*cos(y) - sin(x)*sin(y)
    assert cos(x - y).expand(trig=True) == cos(x)*cos(y) + sin(x)*sin(y)
    assert cos(y - x).expand(trig=True) == cos(x)*cos(y) + sin(x)*sin(y)
    assert cos(2*x).expand(trig=True) == 2*cos(x)**2 - 1
    assert cos(3*x).expand(trig=True) == 4*cos(x)**3 - 3*cos(x)
    assert cos(4*x).expand(trig=True) == 8*cos(x)**4 - 8*cos(x)**2 + 1
    assert cos(2).expand(trig=True) == 2*cos(1)**2 - 1
    assert cos(3).expand(trig=True) == 4*cos(1)**3 - 3*cos(1)

def test_tan():
    x, y = symbols('x,y')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert tan(nan) == nan

    assert tan(oo*I) == I
    assert tan(-oo*I) == -I

    assert tan(0) == 0

    assert tan(atan(x)) == x
    assert tan(asin(x)) == x / sqrt(1 - x**2)
    assert tan(acos(x)) == sqrt(1 - x**2) / x
    assert tan(acot(x)) == 1 / x
    assert tan(atan2(y, x)) == y/x

    assert tan(pi*I) == tanh(pi)*I
    assert tan(-pi*I) == -tanh(pi)*I
    assert tan(-2*I) == -tanh(2)*I

    assert tan(pi) == 0
    assert tan(-pi) == 0
    assert tan(2*pi) == 0
    assert tan(-2*pi) == 0
    assert tan(-3*10**73*pi) == 0

    assert tan(pi/2) == zoo
    assert tan(3*pi/2) == zoo

    assert tan(pi/3) == sqrt(3)
    assert tan(-2*pi/3) == sqrt(3)

    assert tan(pi/4) == S.One
    assert tan(-pi/4) == -S.One
    assert tan(17*pi/4) == S.One
    assert tan(-3*pi/4) == S.One

    assert tan(pi/6) == 1/sqrt(3)
    assert tan(-pi/6) == -1/sqrt(3)
    assert tan(7*pi/6) == 1/sqrt(3)
    assert tan(-5*pi/6) == 1/sqrt(3)

    assert tan(x*I) == tanh(x)*I

    assert tan(k*pi) == 0
    assert tan(17*k*pi) == 0

    assert tan(k*pi*I) == tanh(k*pi)*I

    assert tan(r).is_real == True

    assert tan(10*pi/7) == tan(3*pi/7)
    assert tan(11*pi/7) == -tan(3*pi/7)
    assert tan(-11*pi/7) == tan(3*pi/7)

def test_tan_series():
    x = Symbol('x')
    assert tan(x).series(x, 0, 9) == \
                    x + x**3/3 + 2*x**5/15 + 17*x**7/315 + O(x**9)

def test_tan_rewrite():
    x = Symbol('x')
    neg_exp, pos_exp = exp(-x*I), exp(x*I)
    assert tan(x).rewrite(exp) == I*(neg_exp-pos_exp)/(neg_exp+pos_exp)
    assert tan(x).rewrite(sin) == 2*sin(x)**2/sin(2*x)
    assert tan(x).rewrite(cos) == -cos(x + S.Pi/2)/cos(x)
    assert tan(x).rewrite(cot) == 1/cot(x)
    assert tan(sinh(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, sinh(3)).n()
    assert tan(cosh(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, cosh(3)).n()
    assert tan(tanh(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, tanh(3)).n()
    assert tan(coth(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, coth(3)).n()
    assert tan(sin(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, sin(3)).n()
    assert tan(cos(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, cos(3)).n()
    assert tan(tan(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, tan(3)).n()
    assert tan(cot(x)).rewrite(exp).subs(x, 3).n() == tan(x).rewrite(exp).subs(x, cot(3)).n()
    assert tan(log(x)).rewrite(Pow) == I*(x**-I - x**I)/(x**-I + x**I)

def test_tan_subs():
    x,y = symbols('x,y')
    assert tan(x).subs(tan(x), y) == y
    assert tan(x).subs(x, y) == tan(y)
    assert tan(x).subs(x, S.Pi/2) == zoo
    assert tan(x).subs(x, 3*S.Pi/2) == zoo

def test_cot():
    x, y = symbols('x,y')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert cot(nan) == nan

    assert cot(oo*I) == -I
    assert cot(-oo*I) == I

    assert cot(0) == zoo
    assert cot(2*pi) == zoo

    assert cot(acot(x)) == x
    assert cot(atan(x)) == 1 / x
    assert cot(asin(x)) == sqrt(1 - x**2) / x
    assert cot(acos(x)) == x / sqrt(1 - x**2)
    assert cot(atan2(y,x)) == x/y

    assert cot(pi*I) == -coth(pi)*I
    assert cot(-pi*I) == coth(pi)*I
    assert cot(-2*I) == coth(2)*I

    assert cot(pi) == cot(2*pi) == cot(3*pi)
    assert cot(-pi) == cot(-2*pi) == cot(-3*pi)

    assert cot(pi/2) == 0
    assert cot(-pi/2) == 0
    assert cot(5*pi/2) == 0
    assert cot(7*pi/2) == 0

    assert cot(pi/3) == 1/sqrt(3)
    assert cot(-2*pi/3) == 1/sqrt(3)

    assert cot(pi/4) == S.One
    assert cot(-pi/4) == -S.One
    assert cot(17*pi/4) == S.One
    assert cot(-3*pi/4) == S.One

    assert cot(pi/6) == sqrt(3)
    assert cot(-pi/6) == -sqrt(3)
    assert cot(7*pi/6) == sqrt(3)
    assert cot(-5*pi/6) == sqrt(3)

    assert cot(x*I) == -coth(x)*I
    assert cot(k*pi*I) == -coth(k*pi)*I

    assert cot(r).is_real == True

    assert cot(10*pi/7) == cot(3*pi/7)
    assert cot(11*pi/7) == -cot(3*pi/7)
    assert cot(-11*pi/7) == cot(3*pi/7)

def test_cot_series():
    x = Symbol('x')
    assert cot(x).series(x, 0, 9) == \
                    1/x - x/3 - x**3/45 - 2*x**5/945 - x**7/4725 + O(x**9)

def test_cot_rewrite():
    x = Symbol('x')
    neg_exp, pos_exp = exp(-x*I), exp(x*I)
    assert cot(x).rewrite(exp) == I*(pos_exp+neg_exp)/(pos_exp-neg_exp)
    assert cot(x).rewrite(sin) == 2*sin(2*x)/sin(x)**2
    assert cot(x).rewrite(cos) == -cos(x)/cos(x + S.Pi/2)
    assert cot(x).rewrite(tan) == 1/tan(x)
    assert cot(sinh(x)).rewrite(exp).subs(x, 3).n() == cot(x).rewrite(exp).subs(x, sinh(3)).n()
    assert cot(cosh(x)).rewrite(exp).subs(x, 3).n() == cot(x).rewrite(exp).subs(x, cosh(3)).n()
    assert cot(tanh(x)).rewrite(exp).subs(x, 3).n() == cot(x).rewrite(exp).subs(x, tanh(3)).n()
    assert cot(coth(x)).rewrite(exp).subs(x, 3).n() == cot(x).rewrite(exp).subs(x, coth(3)).n()
    assert cot(sin(x)).rewrite(exp).subs(x, 3).n() == cot(x).rewrite(exp).subs(x, sin(3)).n()
    assert cot(tan(x)).rewrite(exp).subs(x, 3).n() == cot(x).rewrite(exp).subs(x, tan(3)).n()
    assert cot(log(x)).rewrite(Pow) == -I*(x**-I + x**I)/(x**-I - x**I)

def test_cot_subs():
    x,y = symbols('x,y')
    assert cot(x).subs(cot(x), y) == y
    assert cot(x).subs(x, y) == cot(y)
    assert cot(x).subs(x, 0) == zoo
    assert cot(x).subs(x, S.Pi) == zoo

def test_asin():
    x = Symbol('x')

    assert asin(nan) == nan

    assert asin(oo) == -I*oo
    assert asin(-oo) == I*oo

    # Note: asin(-x) = - asin(x)
    assert asin(0)  == 0
    assert asin(1)  == pi/2
    assert asin(-1)  == -pi/2
    assert asin(sqrt(3)/2) == pi/3
    assert asin(-sqrt(3)/2) == -pi/3
    assert asin(sqrt(2)/2) == pi/4
    assert asin(-sqrt(2)/2) == -pi/4
    assert asin(sqrt((5-sqrt(5))/8)) == pi/5
    assert asin(-sqrt((5-sqrt(5))/8)) == -pi/5
    assert asin(Rational(1,2)) == pi/6
    assert asin(-Rational(1,2)) == -pi/6
    assert asin((sqrt(2-sqrt(2)))/2) == pi/8
    assert asin(-(sqrt(2-sqrt(2)))/2) == -pi/8
    assert asin((sqrt(5)-1)/4) == pi/10
    assert asin(-(sqrt(5)-1)/4) == -pi/10
    assert asin((sqrt(3)-1)/sqrt(2**3)) == pi/12
    assert asin(-(sqrt(3)-1)/sqrt(2**3)) == -pi/12

    assert asin(x).diff(x) ==  1/sqrt(1-x**2)

    assert asin(0.2).is_real == True
    assert asin(-2).is_real == False

    assert asin(-2*I) == -I*asinh(2)

def test_asin_series():
    x = Symbol('x')
    assert asin(x).series(x, 0, 9) == \
                    x + x**3/6 + 3*x**5/40 + 5*x**7/112 + O(x**9)
    t5 = asin(x).taylor_term(5, x)
    assert t5 == 3*x**5/40
    assert asin(x).taylor_term(7, x, t5, 0) == 5*x**7/112

def test_asin_rewrite():
    x = Symbol('x')
    assert asin(x).rewrite(log) == -I*log(I*x + sqrt(1 - x**2))
    assert asin(x).rewrite(atan) == 2*atan(x/(1 + sqrt(1 - x**2)))
    assert asin(x).rewrite(acos) == S.Pi/2 - acos(x)

def test_acos():
    x = Symbol('x')
    r = Symbol('r', real=True)

    assert acos(nan) == nan
    assert acos(oo) == I*oo
    assert acos(-oo) == -I*oo

    # Note: acos(-x) = pi - acos(x)
    assert acos(0)  == pi/2
    assert acos(Rational(1,2)) == pi/3
    assert acos(-Rational(1,2)) ==  (2*pi)/3
    assert acos(1)  == 0
    assert acos(-1) == pi
    assert acos(sqrt(2)/2) == pi/4
    assert acos(-sqrt(2)/2) == (3*pi)/4

    assert acos(x).diff(x) == -1/sqrt(1-x**2)

    assert acos(0.2).is_real == True
    assert acos(-2).is_real == False

def test_acos_series():
    x = Symbol('x')
    assert acos(x).series(x, 0, 8) == \
            pi/2 - x - x**3/6 - 3*x**5/40 - 5*x**7/112 + O(x**8)
    assert acos(x).series(x, 0, 8) == pi/2 - asin(x).series(x, 0, 8)
    t5 = acos(x).taylor_term(5, x)
    assert t5 == -3*x**5/40
    assert acos(x).taylor_term(7, x, t5, 0) == -5*x**7/112

def test_acos_rewrite():
    x = Symbol('x')
    assert acos(x).rewrite(log) == pi/2 + I*log(I*x + sqrt(1 - x**2))
    assert acos(0).rewrite(atan) == S.Pi/2
    assert acos(0.5).rewrite(atan) == acos(0.5).rewrite(log)
    assert acos(x).rewrite(asin) == S.Pi/2 - asin(x)

def test_atan():
    x = Symbol('x')

    r = Symbol('r', real=True)

    assert atan(nan) == nan

    assert atan(oo) == pi/2
    assert atan(-oo) == -pi/2

    assert atan(0)  == 0
    assert atan(1)  == pi/4
    assert atan(sqrt(3)) == pi/3
    assert atan(oo) == pi/2
    assert atan(x).diff(x) ==  1/(1+x**2)

    assert atan(r).is_real == True

    assert atan(-2*I) == -I*atanh(2)

def test_atan_rewrite():
    x = Symbol('x')
    assert atan(x).rewrite(log) == I*log((1 - I*x)/(1 + I*x))/2

def test_atan2():
    assert atan2(0, 0) == S.NaN
    assert atan2(0, 1) == 0
    assert atan2(1, 0) == pi/2
    assert atan2(1, -1) == 3*pi/4
    assert atan2(-1, 1) == -pi/4
    assert atan2(0, -1) == pi

def test_acot():
    x = Symbol('x')

    r = Symbol('r', real=True)

    assert acot(nan) == nan

    assert acot(-oo) == 0
    assert acot(oo) == 0
    assert acot(1)  == pi/4
    assert acot(0)  == pi/2
    assert acot(sqrt(3)/3) == pi/3
    assert acot(1/sqrt(3)) == pi/3
    assert acot(-1/sqrt(3)) == -pi/3
    assert acot(x).diff(x) == -1/(1+x**2)

    assert acot(r).is_real == True

    assert acot(I*pi) == -I*acoth(pi)
    assert acot(-2*I) == I*acoth(2)

def test_acot_rewrite():
    x = Symbol('x')
    assert acot(x).rewrite(log) == I*log((x - I)/(x + I))/2

def test_attributes():
    x = Symbol('x')
    assert sin(x).args == (x,)

def test_sincos_rewrite():
    x = Symbol("x")
    y = Symbol("y")
    assert sin(pi/2-x) == cos(x)
    assert sin(pi-x) == sin(x)
    assert cos(pi/2-x) == sin(x)
    assert cos(pi-x) == -cos(x)

def _check_even_rewrite(func, arg):
    """Checks that the expr has been rewritten using f(-x) -> f(x)
    arg : -x
    """
    return func(arg).args[0] == -arg

def _check_odd_rewrite(func, arg):
    """Checks that the expr has been rewritten using f(-x) -> -f(x)
    arg : -x
    """
    return func(arg).func.is_Mul

def _check_no_rewrite(func, arg):
    """Checks that the expr is not rewritten"""
    return func(arg).args[0] == arg

def test_evenodd_rewrite():
    x, y = symbols('x,y')
    a = cos(2) #negative
    b = sin(1) #positive
    even = [cos]
    odd = [sin, tan, cot, asin, atan, acot]
    with_minus = [-1, -2**1024 * E, -pi/105, -x*y, -x-y]
    for func in even:
        for expr in with_minus:
            assert _check_even_rewrite(func, expr)
        assert _check_no_rewrite(func, a*b)
        assert func(x-y) == func(y-x)   #it doesn't matter which form is canonical
    for func in odd:
        for expr in with_minus:
            assert _check_odd_rewrite(func, expr)
        assert _check_no_rewrite(func, a*b)
        assert func(x-y) == -func(y-x)  #it doesn't matter which form is canonical

def test_issue1448():
    x = Symbol('x')
    assert cot(x).inverse() == acot
    assert sin(x).rewrite(cot) == 2*cot(x/2)/(1 + cot(x/2)**2)
    assert cos(x).rewrite(cot) == -(1 - cot(x/2)**2)/(1 + cot(x/2)**2)
    assert tan(x).rewrite(cot) == 1/cot(x)
    assert cot(x).fdiff() == -1 - cot(x)**2

def test_as_leading_term_issue2173():
    x = Symbol('x')
    assert sin(x).as_leading_term(x) == x
    assert cos(x).as_leading_term(x) == 1
    assert tan(x).as_leading_term(x) == x
    assert cot(x).as_leading_term(x) == 1/x
    assert asin(x).as_leading_term(x) == x
    assert acos(x).as_leading_term(x) == x
    assert atan(x).as_leading_term(x) == x
    assert acot(x).as_leading_term(x) == x

def test_leading_terms():
    x = Symbol('x')
    for func in [sin, cos, tan, cot, asin, acos, atan, acot]:
        for arg in (1/x, S.Half):
            eq = func(arg)
            assert eq.as_leading_term(x) == eq

def test_atan2_expansion():
    x, y = symbols("x,y")
    assert cancel(atan2(x+1,x**2).diff(x) - atan((x+1)/x**2).diff(x)) == 0
    assert cancel(atan(x/y).series(x, 0, 5) - atan2(x, y).series(x, 0, 5) \
                  + atan2(0, y) - atan(0)) == O(x**5)
    assert cancel(atan(x/y).series(y, 1, 4) - atan2(x, y).series(y, 1, 4)  \
                  + atan2(x, 1) - atan(x)) == O(y**4)
    assert cancel(atan((x+y)/y).series(y, 1, 3) - atan2(x+y, y).series(y, 1, 3) \
                  + atan2(1+x, 1) - atan(1+x)) == O(y**3)
    assert Matrix([atan2(x, y)]).jacobian([x, y]) \
                  == Matrix([[y/(x**2+y**2), -x/(x**2+y**2)]])

def test_aseries():
    x = Symbol('x')
    def t(n, v, d, e):
        assert abs(n(1/v).evalf() - n(1/x).series(x, dir=d).removeO().subs(x, v)) < e
    t(atan, 0.1, '+', 1e-5)
    t(atan, -0.1, '-', 1e-5)
    t(acot, 0.1, '+', 1e-5)
    t(acot, -0.1, '-', 1e-5)

def test_issue_1321():
    i = Symbol('i', integer=True)
    e = Symbol('e', even=True)
    o = Symbol('o', odd=True)
    x = Symbol('x')

    # unknown parity for variable
    assert cos(4*i*pi) == 1
    assert sin(4*i*pi) == 0
    assert tan(4*i*pi) == 0
    assert cot(4*i*pi) == zoo

    assert cos(3*i*pi) == cos(pi*i) # +/-1
    assert sin(3*i*pi) == 0
    assert tan(3*i*pi) == 0
    assert cot(3*i*pi) == zoo

    assert cos(4.0*i*pi) == 1
    assert sin(4.0*i*pi) == 0
    assert tan(4.0*i*pi) == 0
    assert cot(4.0*i*pi) == zoo

    assert cos(3.0*i*pi) == cos(pi*i) # +/-1
    assert sin(3.0*i*pi) == 0
    assert tan(3.0*i*pi) == 0
    assert cot(3.0*i*pi) == zoo

    assert cos(4.5*i*pi) == cos(0.5*pi*i)
    assert sin(4.5*i*pi) == sin(0.5*pi*i)
    assert tan(4.5*i*pi) == tan(0.5*pi*i)
    assert cot(4.5*i*pi) == cot(0.5*pi*i)

    # parity of variable is known
    assert cos(4*e*pi) == 1
    assert sin(4*e*pi) == 0
    assert tan(4*e*pi) == 0
    assert cot(4*e*pi) == zoo

    assert cos(3*e*pi) == 1
    assert sin(3*e*pi) == 0
    assert tan(3*e*pi) == 0
    assert cot(3*e*pi) == zoo

    assert cos(4.0*e*pi) == 1
    assert sin(4.0*e*pi) == 0
    assert tan(4.0*e*pi) == 0
    assert cot(4.0*e*pi) == zoo

    assert cos(3.0*e*pi) == 1
    assert sin(3.0*e*pi) == 0
    assert tan(3.0*e*pi) == 0
    assert cot(3.0*e*pi) == zoo

    assert cos(4.5*e*pi) == cos(0.5*pi*e)
    assert sin(4.5*e*pi) == sin(0.5*pi*e)
    assert tan(4.5*e*pi) == tan(0.5*pi*e)
    assert cot(4.5*e*pi) == cot(0.5*pi*e)

    assert cos(4*o*pi) == 1
    assert sin(4*o*pi) == 0
    assert tan(4*o*pi) == 0
    assert cot(4*o*pi) == zoo

    assert cos(3*o*pi) == -1
    assert sin(3*o*pi) == 0
    assert tan(3*o*pi) == 0
    assert cot(3*o*pi) == zoo

    assert cos(4.0*o*pi) == 1
    assert sin(4.0*o*pi) == 0
    assert tan(4.0*o*pi) == 0
    assert cot(4.0*o*pi) == zoo

    assert cos(3.0*o*pi) == -1
    assert sin(3.0*o*pi) == 0
    assert tan(3.0*o*pi) == 0
    assert cot(3.0*o*pi) == zoo

    assert cos(4.5*o*pi) == cos(0.5*pi*o)
    assert sin(4.5*o*pi) == sin(0.5*pi*o)
    assert tan(4.5*o*pi) == tan(0.5*pi*o)
    assert cot(4.5*o*pi) == cot(0.5*pi*o)

    # x could be imaginary
    assert cos(4*x*pi) == cos(4*pi*x)
    assert sin(4*x*pi) == sin(4*pi*x)
    assert tan(4*x*pi) == tan(4*pi*x)
    assert cot(4*x*pi) == cot(4*pi*x)

    assert cos(3*x*pi) == cos(3*pi*x)
    assert sin(3*x*pi) == sin(3*pi*x)
    assert tan(3*x*pi) == tan(3*pi*x)
    assert cot(3*x*pi) == cot(3*pi*x)

    assert cos(4.0*x*pi) == cos(4.0*pi*x)
    assert sin(4.0*x*pi) == sin(4.0*pi*x)
    assert tan(4.0*x*pi) == tan(4.0*pi*x)
    assert cot(4.0*x*pi) == cot(4.0*pi*x)

    assert cos(3.0*x*pi) == cos(3.0*pi*x)
    assert sin(3.0*x*pi) == sin(3.0*pi*x)
    assert tan(3.0*x*pi) == tan(3.0*pi*x)
    assert cot(3.0*x*pi) == cot(3.0*pi*x)

    assert cos(4.5*x*pi) == cos(4.5*pi*x)
    assert sin(4.5*x*pi) == sin(4.5*pi*x)
    assert tan(4.5*x*pi) == tan(4.5*pi*x)
    assert cot(4.5*x*pi) == cot(4.5*pi*x)

def test_inverses():
    x = Symbol('x')
    for pair in [[sin, asin], [cos, acos], [tan, atan], [cot, acot]]:
        assert pair[0](x).inverse() == pair[1]

def test_real_imag():
    a,b = symbols('a,b', real=True)
    z = a+b*I
    for deep in [True, False]:
        assert sin(z).as_real_imag(deep=deep) == (sin(a)*cosh(b), cos(a)*sinh(b))
        assert cos(z).as_real_imag(deep=deep) == (cos(a)*cosh(b), -sin(a)*sinh(b))
        assert tan(z).as_real_imag(deep=deep) == (sin(a)*cos(a)/(cos(a)**2+sinh(b)**2), sinh(b)*cosh(b)/(cos(a)**2+sinh(b)**2))
        assert cot(z).as_real_imag(deep=deep) == (sin(a)*cos(a)/(sin(a)**2+sinh(b)**2), -sinh(b)*cosh(b)/(sin(a)**2+sinh(b)**2))
        assert sin(a).as_real_imag(deep=deep) == (sin(a), 0)
        assert cos(a).as_real_imag(deep=deep) == (cos(a), 0)
        assert tan(a).as_real_imag(deep=deep) == (tan(a), 0)
        assert cot(a).as_real_imag(deep=deep) == (cot(a), 0)

@XFAIL
def test_sin_cos_with_infinity():
    # Test for issue 2097
    # http://code.google.com/p/sympy/issues/detail?id=2097
    assert sin(oo) == S.NaN
    assert cos(oo) == S.NaN
