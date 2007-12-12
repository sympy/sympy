from sympy import Function, dsolve, Symbol, sin, cos, sinh, cosh, I, exp, log, \
        simplify, normal, together, ratsimp, powsimp, fraction, radsimp
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
    eq = f(x).diff(x) == 0
    assert dsolve(eq.lhs, f(x)) == Symbol("C1")
    assert dsolve(eq, f(x)) == Symbol("C1")
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode2():
    f = Function("f")
    eq = 3*f(x).diff(x) - 5 == 0
    assert dsolve(eq.lhs, f(x)) == Symbol("C1")+5*x/3
    assert dsolve(eq, f(x)) == Symbol("C1")+5*x/3
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode3():
    f = Function("f")
    eq = 3*f(x).diff(x) == 5
    assert dsolve(eq, f(x)) == Symbol("C1")+5*x/3
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode4():
    f = Function("f")
    eq = 9*f(x).diff(x).diff(x) + f(x) == 0
    assert dsolve(eq, f(x)) == Symbol("C1")*sin(x/3) + Symbol("C2")*cos(x/3)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode5():
    f = Function("f")
    eq = 9*f(x).diff(x).diff(x) == f(x)
    assert dsolve(eq, f(x)) == I*Symbol("C1")*sinh(x/3) + Symbol("C2")*cosh(x/3)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode6():
    f = Function("f")
    # type: (x*exp(-f(x)))'' == 0
    eq = (x*exp(-f(x))).diff(x).diff(x) == 0
    assert dsolve(eq, f(x)) == -log(Symbol("C1")+Symbol("C2")/x)
    checksol(eq, f(x), dsolve(eq, f(x)))

def test_ode7():
    f = Function("f")
    # type: (x*exp(f(x)))'' == 0
    eq = (x*exp(f(x))).diff(x).diff(x) == 0
    assert dsolve(eq, f(x)) == log(Symbol("C1")+Symbol("C2")/x)
    checksol(eq, f(x), dsolve(eq, f(x)))
