from sympy import Function, dsolve, Symbol, sin, cos, sinh, acos, tan, cosh, \
        I, exp, log, simplify, normal, together, ratsimp, powsimp, \
        fraction, radsimp, Eq, sqrt, pi, erf, diff, Rational, asinh, trigsimp
from sympy.abc import x, y, z
from sympy.solvers import deriv_degree
from sympy.solvers.solvers import homogeneous_order
from sympy.utilities.pytest import XFAIL, skip

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
        return simplify(diff(sol.lhs,x,n)-diff(sol.rhs,x,n) - eq) == 0\
        or simplify(trigsimp(diff(sol.lhs,x,n)-diff(sol.rhs,x,n)) - trigsimp(eq)) == 0
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
    eq1 = sin(x)*cos(f(x))+cos(x)*sin(f(x))*f(x).diff(x)
    eq2 = (2*x*f(x)+1)/f(x)+(f(x)-x)/f(x)**2*f(x).diff(x)
    eq3 = 2*x+f(x)*cos(x)+(2*f(x)+sin(x)-sin(f(x)))*f(x).diff(x)
    sol1 = dsolve(eq1,f(x))
    sol2 = dsolve(eq2,f(x))
    sol3 = dsolve(eq3,f(x))
    assert sol1 == Eq(f(x),acos((-C1)/cos(x)))
    assert sol2 == Eq(log(f(x))+x/f(x)+x**2,C1)
    assert sol3 == Eq(f(x)*sin(x)+cos(f(x))+x**2+f(x)**2,C1)
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3)
@XFAIL
def test_ode16():
    # This relates to Issue 1425.  The error is in solve, not dsolve.
    eq = cos(f(x))-(x*sin(f(x))-f(x)**2)*f(x).diff(x)
    sol = dsolve(eq1,f(x))
    assert sol == Eq(x*cos(f(x))+f(x)**3/3,C1)
    assert checksol(eq, f(x), sol)

@XFAIL
def test_ode_exact():
    """
    This is an exact equation that fails under the exact engine. It is caught
    by first order homogeneous albiet with a much contorted solution.  The
    exact engine fails because of a poorly simplified integral of q(0,y)dy,
    where q is the function multiplying f'.  The solutions should be
    Eq((x**2+f(x)**2)**Rational(3,2)+y**3, C1).  The equation below is
    equivalent, but it is so complex that checksol fails, and takes a long
    to do so.
    """
    skip("takes too much time")
    eq = x*sqrt(x**2+f(x)**2)-(x**2*f(x)/(f(x)-sqrt(x**2+f(x)**2)))*f(x).diff(x)
    sol = dsolve(eq, f(x))
    assert sol == Eq(log(x),C1 - 9*(1 + f(x)**2/x**2)**Rational(1,2)*asinh(f(x)/x)/(-27*f(x)/x + 27*(1 + f(x)**2/x**2)**Rational(1,2)) - 9*(1 + f(x)**2/x**2)**Rational(1,2)*log(1 - (1 + f(x)**2/x**2)**Rational(1,2)*f(x)/x + 2*f(x)**2/x**2)/(-27*f(x)/x + 27*(1 + f(x)**2/x**2)**Rational(1,2)) + 9*asinh(f(x)/x)*f(x)/(x*(-27*f(x)/x + 27*(1 + f(x)**2/x**2)**Rational(1,2))) + 9*f(x)*log(1 - (1 + f(x)**2/x**2)**Rational(1,2)*f(x)/x + 2*f(x)**2/x**2)/(x*(-27*f(x)/x + 27*(1 + f(x)**2/x**2)**Rational(1,2))))
    assert checksol(eq, f(x), sol)

def test_homogeneous_order():
    assert homogeneous_order(exp(y/x)+tan(y/x), x, y) == 0
    assert homogeneous_order(x**2 + sin(x)*cos(y), x, y) == None
    assert homogeneous_order(x-y-x*sin(y/x), x, y) == 1
    assert homogeneous_order((x*y+sqrt(x**4+y**4)+x**2*(log(x)-log(y)))/(pi*x**Rational(2,3)*y**Rational(3,2)), x, y) == Rational(-1,6)
    assert homogeneous_order(y/x*cos(y/x)-x/y*sin(y/x)+cos(y/x), x, y) == 0
    assert homogeneous_order(f(x), x, f(x)) == 1
    assert homogeneous_order(f(x)**2, x, f(x)) == 2
    assert homogeneous_order(x*y*z, x, y) == 2
    assert homogeneous_order(x*y*z, x, y, z) == 3
    assert homogeneous_order(x**2*f(x)/sqrt(x**2+f(x)**2), f(x)) == None
    assert homogeneous_order(f(x,y)**2, x, f(x,y), y) == 2
    assert homogeneous_order(f(x,y)**2, x, f(x), y) == None
    assert homogeneous_order(f(x,y)**2, x, f(x,y)) == None
    assert homogeneous_order(f(y,x)**2, x, y, f(x, y)) == None
    assert homogeneous_order(f(y), f(x), x) == None
    assert homogeneous_order(-f(x)/x+1/sin(f(x)/ x), f(x), x) == 0

def test_homogeneous_order_ode1():
    # Type: First order homogeneous, y'=f(y/x)
    eq1 = f(x)/x*cos(f(x)/x)-(x/f(x)*sin(f(x)/x)+cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x)-f(x)-x*sin(f(x)/x)
    eq3 = f(x)+(x*log(f(x)/x)-2*x)*diff(f(x),x)
    eq4 = 2*f(x)*exp(x/f(x))+f(x)*f(x).diff(x)-2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x)+f(x)**3+(x*f(x)**2-2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x)-f(x)*sin(f(x)/x)+x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x+sqrt(f(x)**2-x*f(x)))*f(x).diff(x)-f(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    sol3 = Eq(-f(x)/(1+log(x/f(x))),C1)
    sol4 = Eq(log(f(x)/C1) + 2*exp(x/f(x)), 0)
    sol5 = Eq(log(x*2**Rational(1,2)*(1/x)**Rational(1,2)*f(x)**Rational(1,2)/C1) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(x/C1) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol7 = Eq(log(f(x)/C1) + 2*(1 - x/f(x))**Rational(1,2), 0)
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert dsolve(eq4, f(x)) == sol4
    assert dsolve(eq5, f(x)) == sol5
    assert dsolve(eq6, f(x)) == sol6
    assert dsolve(eq7, f(x)) == sol7

@XFAIL
def test_homogeneous_order_ode1_sol():
    # These are the checksols from test_homogeneous_order_ode1, which fail.
    eq1 = f(x)/x*cos(f(x)/x)-(x/f(x)*sin(f(x)/x)+cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x)-f(x)-x*sin(f(x)/x)
    eq3 = f(x)+(x*log(f(x)/x)-2*x)*diff(f(x),x)
    eq4 = 2*f(x)*exp(x/f(x))+f(x)*f(x).diff(x)-2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x)+f(x)**3+(x*f(x)**2-2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x)-f(x)*sin(f(x)/x)+x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x+sqrt(f(x)**2-x*f(x)))*f(x).diff(x)-f(x)
    eq8 = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    sol3 = Eq(-f(x)/(1+log(x/f(x))),C1)
    sol4 = Eq(log(f(x)/C1) + 2*exp(x/f(x)), 0)
    sol5 = Eq(log(x*2**Rational(1,2)*(1/x)**Rational(1,2)*f(x)**Rational(1,2)/C1) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(x/C1) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol7 = Eq(log(f(x)/C1) + 2*(1 - x/f(x))**Rational(1,2), 0)
    sol8 =  Eq(-Integral(-1/((1 - u2**2)**Rational(1,2)*u2), (u2, _a, x/f(x))) + log(f(x)/C1), 0)
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3)
    assert checksol(eq4, f(x), sol4)
    assert checksol(eq5, f(x), sol5)
    assert checksol(eq6, f(x), sol6)
    assert checksol(eq7, f(x), sol7)
    assert checksol(eq8, f(x), sol8)

def test_homogeneous_order_ode2():
    eq1 = f(x).diff(x)-f(x)/x+1/sin(f(x)/x)
    eq2 = x**2+f(x)**2-2*x*f(x)*f(x).diff(x)
    eq3 = x*exp(f(x)/x)+f(x)-x*f(x).diff(x)
    sol1 = Eq(f(x), x*acos(log(x/C1)))
    sol2 = [Eq(f(x), (-C1*x + x**2)**Rational(1,2)), Eq(f(x), -(-C1*x + x**2)**Rational(1,2))]
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2[0])
    assert checksol(eq2, f(x), sol2[1])
    assert checksol(eq3, f(x), sol3)

@XFAIL
def test_homogeneous_order_ode3():
    # This fails because of issue 1429
    eq = x+f(x)-(x-f(x))*f(x).diff(x)
    sol = Eq(f(x), 2*x*atan(C1*x))
    assert dsolve(eq, f(x)) == sol
    assert checksol(eq, f(x), sol)

def test_homogeneous_order_ode4():
    # This can be solved explicitly, but the the integration engine cannot handle
    # it (see issue 1452).  The explicit solution is included in an XFAIL test
    # below. checksol fails for this equation, so its test is in
    # test_homogeneous_order_ode1_sol above. It has to compare string
    # expressions because u2 and _a are dummy variables.
    eq = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    solstr =  "-Integral(-1/((1 - _u2**2)**(1/2)*_u2), (_u2, __a, x/f(x))) + log(f(x)/C1) == 0"
    assert str(dsolve(eq, f(x))) == solstr


@XFAIL
def test_homogeneous_order_ode4_explicit():
    eq = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    sol = Eq(f(x)**2-C1*x, f(x)*sqrt(f(x)**2-x**2))
    assert dsolve(eq, f(x)) == sol
    assert checksol(eq, f(x)) == sol
