from sympy import sin, exp, Function, Symbol, S, Pow, Eq, I, sinh, cosh, acos,\
cos, log, Rational, sqrt, Integral
from sympy.solvers.ode import constantsimp
from sympy.utilities.pytest import XFAIL


x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
C1 = Symbol('C1')
C2 = Symbol('C2')
C3 = Symbol('C3')
f = Function('f')

def test_constant_mul():
    # We want C1 (Constant) below to absorb the y's, but not the x's
    assert constantsimp(y*C1, x, 1) == C1
    assert constantsimp(x*C1, x, 1) == x*C1
    assert constantsimp(C1*y, x, 1) == C1
    assert constantsimp(C1*x, x, 1) == x*C1
    assert constantsimp(2*C1, x, 1) == C1
    assert constantsimp(C1*2, x, 1) == C1
    assert constantsimp(y*C1*x, x, 1) == C1*x
    assert constantsimp(x*y*C1, x, 1) == x*C1
    assert constantsimp(y*x*C1, x, 1) == x*C1
    assert constantsimp(C1*y*(y + 1), x, 1) == C1
    assert constantsimp(y*C1*(y + 1), x, 1) == C1
    assert constantsimp(x*(y*C1), x, 1) == x*C1
    assert constantsimp(x*(C1*y), x, 1) == x*C1
    assert constantsimp(C1*(x*y), x, 1) == C1*x
    assert constantsimp((x*y)*C1, x, 1) == x*C1
    assert constantsimp((y*x)*C1, x, 1) == x*C1
    assert constantsimp(y*(y + 1)*C1, x, 1) == C1
    assert constantsimp(C1*x*y, x, 1) == C1*x
    assert constantsimp(x*C1*y, x, 1) == x*C1
    assert constantsimp((C1*x)*y, x, 1) == C1*x
    assert constantsimp(y*(x*C1), x, 1) == x*C1
    assert constantsimp((x*C1)*y, x, 1) == x*C1
    assert constantsimp(C1*x*y*x*y*2, x, 1) == C1*x**2
    assert constantsimp(C1*x*y*z, x, 1) == C1*x
    assert constantsimp(C1*x*y**2*sin(z), x, 1) == C1*x
    assert constantsimp(C1*C1, x, 1) == C1
    assert constantsimp(C1*C2, x, 2) == C1
    assert constantsimp(C2*C2, x, 2) == C1
    assert constantsimp(C1*C1*C2, x, 2) == C1
    assert constantsimp(C1*x*2**x, x, 1) == C1*x*2**x

def test_constant_add():
    assert constantsimp(C1 + C1, x, 1) == C1
    assert constantsimp(C1 + 2, x, 1) == C1
    assert constantsimp(2 + C1, x, 1) == C1
    assert constantsimp(C1 + y, x, 1) == C1
    assert constantsimp(C1 + x, x, 1) == C1 + x
    assert constantsimp(C1 + x + y + x*y + 2, x, 1) == C1 + x + x*y
    assert constantsimp(C1 + x + 2**x + y + 2, x, 1) == C1 + x + 2**x
    assert constantsimp(C1 + C1, x, 1) == C1
    assert constantsimp(C1 + C2, x, 2) == C1
    assert constantsimp(C2 + C1, x, 2) == C1
    assert constantsimp(C1 + C2 + C1, x, 2) == C1

def test_constant_power_as_base():
    assert constantsimp(C1**C1, x, 1) == C1
    assert constantsimp(Pow(C1,C1), x, 1) == C1
    assert constantsimp(C1**C1, x, 1) == C1
    assert constantsimp(C1**C2, x, 2) == C1
    assert constantsimp(C2**C1, x, 2) == C1
    assert constantsimp(C2**C2, x, 2) == C1
    assert constantsimp(C1**y, x, 1) == C1
    assert constantsimp(C1**x, x, 1) == C1**x
    assert constantsimp(C1**2, x, 1) == C1
    assert constantsimp(C1**(x*y), x, 1) == C1**(x*y)

def test_constant_power_as_exp():
    assert constantsimp(x**C1, x, 1) == x**C1
    assert constantsimp(y**C1, x, 1) == C1
    assert constantsimp(x**y**C1, x, 1) == x**C1
    assert constantsimp((x**y)**C1, x, 1) == (x**y)**C1
    assert constantsimp(x**(y**C1), x, 1) == x**C1
    assert constantsimp(x**C1**y, x, 1) == x**C1
    assert constantsimp(x**(C1**y), x, 1) == x**C1
    assert constantsimp((x**C1)**y, x, 1) == (x**C1)**y
    assert constantsimp(2**C1, x, 1) == C1
    assert constantsimp(S(2)**C1, x, 1) == C1
    assert constantsimp(exp(C1), x, 1) == C1
    assert constantsimp(exp(C1+x), x, 1) == exp(C1+x)
    assert constantsimp(Pow(2, C1), x, 1) == C1

def test_constant_function():
    assert constantsimp(sin(C1), x, 1) == C1
    assert constantsimp(f(C1), x, 1) == C1
    assert constantsimp(f(C1, C1), x, 1) == C1
    assert constantsimp(f(C1, C2), x, 2) == C1
    assert constantsimp(f(C2, C1), x, 2) == C1
    assert constantsimp(f(C2, C2), x, 2) == C1
    assert constantsimp(f(C1, x), x, 1) == f(C1, x)
    assert constantsimp(f(C1, y), x, 1) == C1
    assert constantsimp(f(y, C1), x, 1) == C1
    assert constantsimp(f(C1, y, C2), x, 2) == C1

@XFAIL
def test_constant_function_multiple():
    # The rules to not renumber in this case would be too complicated, and
    # dsolve is not likely to ever encounter anything remotely like this.
    assert constantsimp(f(C1, C1, x), x, 1) == f(C1, C1, x)

def test_constant_multiple():
    assert constantsimp(C1*2 + 2, x, 1) == C1
    assert constantsimp(x*2/C1, x, 1) == C1*x
    assert constantsimp(C1**2*2 + 2, x, 1) == C1
    assert constantsimp(sin(2*C1) + x + sqrt(2), x, 1) == C1 + x
    assert constantsimp(2*C1 + C2, x, 2) == C1

def test_ode_solutions():
    # only a few examples here, the rest will be tested in the actual dsolve tests
    assert constantsimp(C1*exp(2*x)+exp(x)*(C2+C3), x, 3) == C1*exp(x)+C2*exp(2*x)
    assert constantsimp(Eq(f(x),I*C1*sinh(x/3) + C2*cosh(x/3)), x, 2) == \
    Eq(f(x), C1*sinh(x/3) + C2*cosh(x/3))
    assert constantsimp(Eq(f(x),acos((-C1)/cos(x))), x, 1) == \
    Eq(f(x),acos(C1/cos(x)))
    assert constantsimp(Eq(log(f(x)/C1) + 2*exp(x/f(x)), 0), x, 1) == \
    Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    assert constantsimp(Eq(log(x*2**Rational(1,2)*(1/x)**Rational(1,2)*f(x)\
    **Rational(1,2)/C1) + x**2/(2*f(x)**2), 0), x, 1) == \
    Eq(log(C1*x*(1/x)**Rational(1,2)*f(x)**Rational(1,2)) + x**2/(2*f(x)**2), 0)
    assert constantsimp(Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(x/C1) - \
    cos(f(x)/x)*exp(-f(x)/x)/2, 0), x, 1) == Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + \
    log(C1*x) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    u2 = Symbol('u2')
    _a = Symbol('_a')
    assert constantsimp(Eq(-Integral(-1/((1 - u2**2)**Rational(1,2)*u2), \
    (u2, _a, x/f(x))) + log(f(x)/C1), 0), x, 1) == Eq(-Integral(-1/(u2*(1 - \
    u2**2)**Rational(1,2)), (u2, _a, x/f(x))) + log(C1*f(x)), 0)
    assert map(lambda i: constantsimp(i, x, 1), [Eq(f(x), (-C1*x + x**2)**\
    Rational(1,2)), Eq(f(x), -(-C1*x + x**2)**Rational(1,2))]) == [Eq(f(x), \
    (C1*x + x**2)**Rational(1,2)), Eq(f(x), -(C1*x + x**2)**Rational(1,2))]



