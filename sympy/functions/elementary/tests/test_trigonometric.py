from sympy import symbols, Symbol, nan, oo, I, sinh, sin, acot, pi, atan, \
        acos, Rational, sqrt, asin, acot, cot, coth, E, S, tan, tanh, cos, \
        cosh, atan2, exp, asinh, acoth, atanh, O

def test_sin():
    x, y = symbols('xy')

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


def test_cos():
    x, y = symbols('xy')

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

def test_tan():
    x, y = symbols('xy')

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

    assert tan(pi*I) == tanh(pi)*I
    assert tan(-pi*I) == -tanh(pi)*I
    assert tan(-2*I) == -tanh(2)*I

    assert tan(pi) == 0
    assert tan(-pi) == 0
    assert tan(2*pi) == 0
    assert tan(-2*pi) == 0
    assert tan(-3*10**73*pi) == 0

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

def test_tan_rewrite():
    x = Symbol('x')
    neg_exp, pos_exp = exp(-x*I), exp(x*I)
    assert tan(x).rewrite(exp) == I*(neg_exp-pos_exp)/(neg_exp+pos_exp)
    assert tan(x).rewrite(sin) == 2*sin(x)**2/sin(2*x)
    assert tan(x).rewrite(cos) == -cos(x + S.Pi/2)/cos(x)
    assert tan(x).rewrite(cot) == 1/cot(x)


def test_cot():
    x, y = symbols('xy')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert cot(nan) == nan

    assert cot(oo*I) == -I
    assert cot(-oo*I) == I

    assert cot(acot(x)) == x
    assert cot(atan(x)) == 1 / x
    assert cot(asin(x)) == sqrt(1 - x**2) / x
    assert cot(acos(x)) == x / sqrt(1 - x**2)

    assert cot(pi*I) == -coth(pi)*I
    assert cot(-pi*I) == coth(pi)*I
    assert cot(-2*I) == coth(2)*I

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

def test_cot_rewrite():
    x = Symbol('x')
    neg_exp, pos_exp = exp(-x*I), exp(x*I)
    assert cot(x).rewrite(exp) == I*(pos_exp+neg_exp)/(pos_exp-neg_exp)
    assert cot(x).rewrite(sin) == 2*sin(2*x)/sin(x)**2
    assert cot(x).rewrite(cos) == -cos(x)/cos(x + S.Pi/2)
    assert cot(x).rewrite(tan) == 1/tan(x)

def test_asin():
    x = Symbol('x')

    assert asin(nan) == nan

    assert asin(oo) == -I*oo
    assert asin(-oo) == I*oo

    assert asin(0)  == 0
    assert asin(Rational(1,2)) == pi/6
    assert asin(1)  == pi/2
    assert asin(sqrt(3)/2) == pi/3

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

def test_acos():
    x = Symbol('x')
    r = Symbol('r', real=True)

    assert acos(nan) == nan
    assert acos(oo) == I*oo
    assert acos(-oo) == -I*oo

    assert acos(0)  == pi/2
    assert acos(Rational(1,2)) == pi/3
    assert acos(1)  == 0
    assert acos(-1) == pi
    assert acos(sqrt(2)/2) == pi/4
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
    x, y = symbols('xy')
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

