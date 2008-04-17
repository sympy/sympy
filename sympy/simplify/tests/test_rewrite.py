from sympy import sin, cos, exp, cot, sqrt, I, E, pi, symbols, Function
from sympy.simplify import cancel, trim, apart, together
from sympy.integrals import integrate
from sympy.utilities.pytest import XFAIL

x,y,n = symbols('xyn')

def test_has():
    assert cot(x).has(x)
    assert cot(x).has(cot)
    assert not cot(x).has(sin)
    assert sin(x).has(x)
    assert sin(x).has(sin)
    assert not sin(x).has(cot)

def test_sin_exp_rewrite():
    assert sin(x).rewrite(sin, exp) == -I/2*(exp(I*x)-exp(-I*x))

def test_cancel():
    assert cancel((x**2-2)/(x+sqrt(2))) == x - sqrt(2)
    assert cancel((x**2-y**2)/(x-y), x, y) == x + y
    assert cancel((x**2-y)/(x-y)) == 1/(x - y)*(x**2 - y)
    assert cancel((x**2-y**2)/(x-y), x) == x + y
    assert cancel((x**2-y**2)/(x-y), y) == x + y
    assert cancel((x**2-y**2)/(x-y)) == x + y
    assert cancel((E*x+2)/(x-pi)*(x-1)) in  \
            [-1/(x - pi)*(2 - 2*x + E*x - E*x**2),
             1/(x - pi)*(-2 + 2*x - E*x + E*x**2)]
    assert cancel((x**3-1)/(x-1) < sin(x)) == (1 + x + x**2 < sin(x))
    assert cancel((x**2-1)/(x-1) == (x**2+1)/(x-I), x) == (1 + x == I + x)
    assert cancel((x**2-1)/(x-1) + (x**2+1)/(x-I), x) == 1 + I + 2*x
    assert cancel((x**2-1)/(x-1) + (x**2+1)/(x-I), y) in [
            1/(1 - x)*(1 - x**2) + 1/(x - I)*(1 + x**2),
            -1/(1 - x)*(-1 + x**2) + 1/(x - I)*(1 + x**2)]

def test_trim():
    f = Function('f')

    assert trim((f(x)**2+f(x))/f(x)) == 1 + f(x)
    assert trim((sin(x)**2+sin(x))/sin(x)) == 1 + sin(x)

    assert trim((f(x)+y*f(x))/f(x)) == 1 + y

    expr = integrate(1/(x**3+1), x)

    assert trim(together(expr.diff(x))) == 1/(x**3+1)
    assert cancel(together(expr.diff(x))) == 1/(x**3+1)

    expr = together(expr.subs(x, sin(x)).diff(x))

    assert trim(expr) == cos(x)/(1 + sin(x)**3)

    assert trim((2 * (1/n - cos(n * pi)/n))/pi) == \
        1/pi/n*(2 - 2*cos(pi*n))

    assert trim(sin((f(x)**2+f(x))/f(x))) == sin(1 + f(x))

    assert trim(exp(x)*sin(x)/2 + cos(x)*exp(x)) == \
        exp(x)*(sin(x) + 2*cos(x))/2

@XFAIL  # because of #666
def test_trim_xfail():
    assert trim((exp(3*x)+exp(2*x*y) + y*exp(x))/exp(x)) == \
        y + exp(2*x) + exp(-x + 2*x*y)


def test_apart1():
    assert apart(1/(x+2)/(x+1), x) == 1/(1 + x) - 1/(2 + x)

def test_apart2():
    assert apart(1/(x+1)/(x+5), x) == -1/(5 + x)/4 + 1/(1 + x)/4
    assert apart((x**2+1)/(x+1) == sin(x), x) == (x - 1 + 2/(x+1) == sin(x))
    assert apart((E*x+2)/(x-pi)*(x-1), x) in [
            2 - E + E*pi + E*x - 1/(x - pi)*(2 - 2*pi + E*pi - E*pi**2),
            2 - E + E*pi + E*x + 1/(x - pi)*(-2 + 2*pi - E*pi + E*pi**2)]
