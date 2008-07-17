from sympy import Function, dsolve, Symbol, sin, cos, sinh, cosh, I, exp, log, \
        simplify, normal, together, ratsimp, powsimp, fraction, radsimp, Eq, \
        sqrt, pi, erf, diff
from sympy.abc import x

def checksol(eq, func, sol):
    """Substitutes sol for func in eq and checks, that the result is 0

    Only works when func is one function, like f(x) and sol just one
    solution like A*sin(x)+B*cos(x).
    """
    s = eq.subs(func, sol)
    if isinstance(s, bool):
        assert s
    elif s:
        return
    else:
        assert s.rhs == 0
        s = simplify(s.lhs)
        assert s == 0


def test_ode1():
    f = Function("f")
    eq = Eq(f(x).diff(x), 0)
    assert dsolve(eq.lhs, f(x)) == Symbol("C1")
    assert dsolve(eq, f(x)) == Symbol("C1")
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode2():
    f = Function("f")
    eq = Eq(3*f(x).diff(x) - 5, 0)
    assert dsolve(eq.lhs, f(x)) == Symbol("C1")+5*x/3
    assert dsolve(eq, f(x)) == Symbol("C1")+5*x/3
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode3():
    f = Function("f")
    eq = Eq(3*f(x).diff(x), 5)
    assert dsolve(eq, f(x)) == Symbol("C1")+5*x/3
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode4():
    f = Function("f")
    eq = Eq(9*f(x).diff(x, x) + f(x), 0)
    assert dsolve(eq, f(x)) == Symbol("C1")*sin(x/3) + Symbol("C2")*cos(x/3)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode5():
    f = Function("f")
    eq = Eq(9*f(x).diff(x, x), f(x))
    assert dsolve(eq, f(x)) == I*Symbol("C1")*sinh(x/3) + Symbol("C2")*cosh(x/3)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode6():
    f = Function("f")
    # type: (x*exp(-f(x)))'' == 0
    eq = Eq((x*exp(-f(x))).diff(x, x), 0)
    assert dsolve(eq, f(x)) == -log(Symbol("C1")+Symbol("C2")/x)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode7():
    f = Function("f")
    # type: (x*exp(f(x)))'' == 0
    eq = Eq((x*exp(f(x))).diff(x, x), 0)
    assert dsolve(eq, f(x)) == log(Symbol("C1")+Symbol("C2")/x)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode8():
    f = Function("f")
    # type: a(x)f'(x)+b(x)*f(x)+c(x)=0
    eq = Eq(x**2*f(x).diff(x) + 3*x*f(x) - sin(x)/x, 0)
    assert dsolve(eq, f(x)) == (Symbol("C1")-cos(x))/x**3
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode9():
    f = Function("f")
    # type:first order linear form f'(x)+p(x)f(x)=q(x)
    eq = Eq(f(x).diff(x) + x*f(x), x**2)
    assert dsolve(eq, f(x)) == exp(-x**2/2)*(sqrt(2)*sqrt(pi)*I*erf(I*x/sqrt(2))/2 + x*exp(x**2/2) + Symbol("C1"))
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode10():
    f = Function("f")
    #type:2nd order, constant coefficients (two real different roots)
    eq = Eq(f(x).diff(x,x) - 3*diff(f(x),x) + 2*f(x), 0)
    assert dsolve(eq, f(x)) in [
        Symbol("C1")*exp(2*x) + Symbol("C2")*exp(x),
        Symbol("C1")*exp(x) + Symbol("C2")*exp(2*x),
    ]
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode11():
    f = Function("f")
    #type:2nd order, constant coefficients (two real equal roots)
    eq = Eq(f(x).diff(x,x) - 4*diff(f(x),x) + 4*f(x), 0)
    assert dsolve(eq, f(x)) == (Symbol("C1") + Symbol("C2")*x)*exp(2*x)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode12():
    f = Function("f")
    #type:2nd order, constant coefficients (two complex roots)
    eq = Eq(f(x).diff(x,x)+2*diff(f(x),x)+3*f(x), 0)
    assert dsolve(eq, f(x)) == (Symbol("C1")*sin(x*sqrt(2))+Symbol("C2")*cos(x*sqrt(2)))*exp(-x)
    checksol(eq, f(x), dsolve(eq, f(x)))

