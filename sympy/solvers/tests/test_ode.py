from sympy import Function, dsolve, Symbol, sin, cos, sinh, cosh, I, exp, log, \
        simplify, normal, together, ratsimp, powsimp, fraction, radsimp, Eq, \
        sqrt, pi, erf, diff, acos
from sympy.abc import x
from sympy.solvers import deriv_degree

C1 = Symbol('C1')
C2 = Symbol('C2')
f = Function('f')

# Note that if the ODE solver changes, these tests could fail but still be
# correct because the arbitrary constants could represent different constants.
# See issue 1336.

def checksol(eq, func, sol):
    """Substitutes sol for func in eq and checks that the result is 0.

    Only works when func is one function, like f(x) and sol just one
    solution like A*sin(x)+B*cos(x).

    It attempts to substitute the solution for f in the original equation.
    If it can't do that, it takes n derivatives of the solution, where eq is of
    order n and checks to see if that is equal to the solution.

    Returns True if the solution checks and False otherwise.
    """

    if sol.lhs == func:
            s = eq.subs(func,sol.rhs)
    elif sol.rhs == func:
            s = eq.subs(func,sol.lhs)
    else:
        # If we cannot substitute f, try seeing if the nth derivative is equal
        n = deriv_degree(eq, func)
        return simplify(diff(sol.lhs,x,n)-diff(sol.rhs,x,n)) == \
        simplify(eq)
    if s == 0:
        return True
    if isinstance(s, bool):
        return s
    elif s:
        return s
    else:
        test1 = (s.rhs == 0)
        s = simplify(s.lhs)
        test2 = (s == 0)
        return test1 and test2


def test_ode1():
    eq = Eq(f(x).diff(x), 0)
    sol1 = dsolve(eq.lhs, f(x))
    sol2 = dsolve(eq, f(x))
    assert sol1 == Eq(f(x),C1)
    assert sol2 == Eq(f(x),C1)
    assert checksol(eq, f(x), sol2)

def test_ode2():
    eq = Eq(3*f(x).diff(x) - 5, 0)
    sol1 = dsolve(eq.lhs, f(x))
    sol2 = dsolve(eq, f(x))
    assert sol1 == Eq(f(x),C1+5*x/3)
    assert sol2 == Eq(f(x),C1+5*x/3)
    assert checksol(eq, f(x), sol2)

def test_ode3():
    eq = Eq(3*f(x).diff(x), 5)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),C1+5*x/3)
    assert checksol(eq, f(x), sol)

def test_ode4():
    eq = Eq(9*f(x).diff(x, x) + f(x), 0)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),C1*sin(x/3) + C2*cos(x/3))
    assert checksol(eq, f(x), sol)

def test_ode5():
    eq = Eq(9*f(x).diff(x, x), f(x))
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),I*C1*sinh(x/3) + C2*cosh(x/3))
    assert checksol(eq, f(x), sol)

def test_ode6():
    # Type: (x*exp(-f(x)))'' == 0
    eq = Eq((x*exp(-f(x))).diff(x, x), 0)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),-log(C1+C2/x))
    assert checksol(eq, f(x), sol)

def test_ode7():
    # Type: (x*exp(f(x)))'' == 0
    eq = Eq((x*exp(f(x))).diff(x, x), 0)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),log(C1+C2/x))
    assert checksol(eq, f(x), sol)

def test_ode8():
    # Type: a(x)f'(x)+b(x)*f(x)+c(x)=0
    eq = Eq(x**2*f(x).diff(x) + 3*x*f(x) - sin(x)/x, 0)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),(C1-cos(x))/x**3)
    assert checksol(eq, f(x), sol)

def test_ode9():
    # Type: first order linear form f'(x)+p(x)f(x)=q(x)
    eq = Eq(f(x).diff(x) + x*f(x), x**2)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),exp(-x**2/2)*(sqrt(2)*sqrt(pi)*I*erf(I*x/sqrt(2))/2 \
    + x*exp(x**2/2) + C1))
    assert checksol(eq, f(x), sol)

def test_ode10():
    # Type: 2nd order, constant coefficients (two real different roots)
    eq = Eq(f(x).diff(x,x) - 3*diff(f(x),x) + 2*f(x), 0)
    sol = dsolve(eq, f(x))
    assert sol in [
        Eq(f(x),C1*exp(2*x) + C2*exp(x)),
        Eq(f(x),C1*exp(x) + C2*exp(2*x)),
    ]
    assert checksol(eq, f(x), sol)

def test_ode11():
    # Type: 2nd order, constant coefficients (two real equal roots)
    eq = Eq(f(x).diff(x,x) - 4*diff(f(x),x) + 4*f(x), 0)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),(C1 + C2*x)*exp(2*x))
    assert checksol(eq, f(x), sol)

def test_ode12():
    # Type: 2nd order, constant coefficients (two complex roots)
    eq = Eq(f(x).diff(x,x)+2*diff(f(x),x)+3*f(x), 0)
    sol = dsolve(eq, f(x))
    assert sol == Eq(f(x),(C1*sin(x*sqrt(2))+C2*cos(x*sqrt(2)))*exp(-x))
    assert checksol(eq, f(x), sol)

def test_ode13():
    eq1 = Eq(3*f(x).diff(x) -1,0)
    eq2 = Eq(x*f(x).diff(x) -1,0)
    sol1 = dsolve(eq1, f(x))
    sol2 = dsolve(eq2, f(x))
    assert sol1 == Eq(f(x),x/3 + C1)
    assert sol2 == Eq(f(x),log(x) + C1)
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)

def test_ode14():
    # Type: Bernoulli, f'(x)+p(x)f(x)=q(x)f(x)**n
    eq = Eq(x*f(x).diff(x)+f(x)-f(x)**2,0)
    sol = dsolve(eq,f(x))
    assert sol == Eq(f(x),1/(x*(C1 + 1/x)))
    assert checksol(eq, f(x), sol)

def test_ode15():
    # Type: Exact differential equation, p(x,f)+q(x,f)f'=0,
    # where dp/dy == dq/dx
    eq1 = cos(f(x))-(x*sin(f(x))-f(x)**2)*f(x).diff(x)
    eq2 = sin(x)*cos(f(x))+cos(x)*sin(f(x))*f(x).diff(x)
    eq3 = (2*x*f(x)+1)/f(x)+(f(x)-x)/f(x)**2*f(x).diff(x)
    eq4 = 2*x+f(x)*cos(x)+(2*f(x)+sin(x)-sin(f(x)))*f(x).diff(x)
    sol1 =  dsolve(eq1,f(x))
    sol2 = dsolve(eq2,f(x))
    sol3 = dsolve(eq3,f(x))
    sol4 = dsolve(eq4,f(x))
    assert sol1 == Eq(x*cos(f(x))+f(x)**3/3,C1)
    assert sol2 == Eq(f(x),acos((-C1)/cos(x)))
    assert sol3 == Eq(log(f(x))+x/f(x)+x**2,C1)
    assert sol4 == Eq(f(x)*sin(x)+cos(f(x))+x**2+f(x)**2,C1)
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3)
