from sympy import Function, dsolve, Symbol, sin, cos, sinh, acos, tan, cosh, \
        I, exp, log, simplify, normal, together, ratsimp, powsimp, \
        fraction, radsimp, Eq, sqrt, pi, erf, diff, Rational, asinh, trigsimp
from sympy.abc import x, y, z
from sympy.solvers import deriv_degree
from sympy.solvers.solvers import homogeneous_order
from sympy.utilities.pytest import XFAIL, skip

C1 = Symbol('C1')
C2 = Symbol('C2')
C3 = Symbol('C3')
C4 = Symbol('C4')
C5 = Symbol('C5')
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
        return simplify(diff(sol.lhs,x,n)-diff(sol.rhs,x,n) - eq) == 0 \
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
    assert sol == Eq(f(x),C1*exp(-x/3) + C2*exp(x/3))
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
    assert sol == Eq(f(x), (C1*sin(x*sqrt(2)) + C2*cos(x*sqrt(2)))*exp(-x))
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
    # Type: Bernoulli, f'(x) + p(x)*f(x) == q(x)*f(x)**n
    eq = Eq(x*f(x).diff(x) + f(x) - f(x)**2,0)
    sol = dsolve(eq,f(x))
    assert sol == Eq(f(x),1/(x*(C1 + 1/x)))
    assert checksol(eq, f(x), sol)

def test_ode15():
    # Type: Exact differential equation, p(x,f) + q(x,f)*f' == 0,
    # where dp/df == dq/dx
    eq1 = sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x)
    eq2 = (2*x*f(x) + 1)/f(x) + (f(x) - x)/f(x)**2*f(x).diff(x)
    eq3 = 2*x + f(x)*cos(x) + (2*f(x) + sin(x) - sin(f(x)))*f(x).diff(x)
    sol1 = dsolve(eq1,f(x))
    sol2 = dsolve(eq2,f(x))
    sol3 = dsolve(eq3,f(x))
#    assert sol1 == Eq(f(x),acos((C1)/cos(x))) # refator into hints engine (separable)
    assert sol2 == Eq(log(f(x))+x/f(x)+x**2,C1)
    assert sol3 == Eq(f(x)*sin(x)+cos(f(x))+x**2+f(x)**2,C1)
#    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3)

@XFAIL
def test_ode16():
    # This relates to Issue 1425.  The error is in solve, not dsolve.
    eq = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
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
    equivalent, but it is so complex that checksol fails, and takes a long time
    to do so.
    """
    skip("takes too much time")
    eq = x*sqrt(x**2 + f(x)**2) - (x**2*f(x)/(f(x) - sqrt(x**2 + f(x)**2)))*f(x).diff(x)
    sol = dsolve(eq, f(x))
    assert sol == Eq(log(x),C1 - 9*sqrt(1 + f(x)**2/x**2)*asinh(f(x)/x)/(-27*f(x)/x + \
    27*sqrt(1 + f(x)**2/x**2)) - 9*sqrt(1 + f(x)**2/x**2)*log(1 - sqrt(1 + f(x)**2/x**2)*\
    f(x)/x + 2*f(x)**2/x**2)/(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2)) \
    + 9*asinh(f(x)/x)*f(x)/(x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))) \
    + 9*f(x)*log(1 - sqrt(1 + f(x)**2/x**2)*f(x)/x + 2*f(x)**2/x**2)/\
    (x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))))
    assert checksol(eq, f(x), sol)

def test_homogeneous_order():
    assert homogeneous_order(exp(y/x) + tan(y/x), x, y) == 0
    assert homogeneous_order(x**2 + sin(x)*cos(y), x, y) == None
    assert homogeneous_order(x - y - x*sin(y/x), x, y) == 1
    assert homogeneous_order((x*y + sqrt(x**4+y**4) + x**2*(log(x) - log(y)))/\
    (pi*x**Rational(2,3)*y**Rational(3,2)), x, y) == Rational(-1,6)
    assert homogeneous_order(y/x*cos(y/x) - x/y*sin(y/x) + cos(y/x), x, y) == 0
    assert homogeneous_order(f(x), x, f(x)) == 1
    assert homogeneous_order(f(x)**2, x, f(x)) == 2
    assert homogeneous_order(x*y*z, x, y) == 2
    assert homogeneous_order(x*y*z, x, y, z) == 3
    assert homogeneous_order(x**2*f(x)/sqrt(x**2 + f(x)**2), f(x)) == None
    assert homogeneous_order(f(x,y)**2, x, f(x,y), y) == 2
    assert homogeneous_order(f(x,y)**2, x, f(x), y) == None
    assert homogeneous_order(f(x,y)**2, x, f(x,y)) == None
    assert homogeneous_order(f(y,x)**2, x, y, f(x, y)) == None
    assert homogeneous_order(f(y), f(x), x) == None
    assert homogeneous_order(-f(x)/x + 1/sin(f(x)/ x), f(x), x) == 0

def test_homogeneous_coeff_ode1():
    # Type: First order homogeneous, y'=f(y/x)
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x),x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    sol3 = Eq(-f(x)/(1+log(x/f(x))),C1)
    sol4 = Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    sol5 = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(C1*x) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol7 = Eq(log(C1*f(x)) + 2*sqrt(1 - x/f(x)), 0)
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert dsolve(eq4, f(x)) == sol4
    assert dsolve(eq5, f(x)) == sol5
#    assert dsolve(eq6, f(x)) == sol6
    assert dsolve(eq7, f(x)) == sol7

@XFAIL
def test_homogeneous_coeff_ode1_sol():
    # These are the checksols from test_homogeneous_coeff_ode1, which fail.
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x),x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    # test_homogeneous_coeff_ode4
    eq8 = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    sol3 = Eq(-f(x)/(1+log(x/f(x))),C1)
    sol4 = Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    sol5 = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(C1*x) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol7 = Eq(log(C1*f(x)) + 2*sqrt(1 - x/f(x)), 0)
    sol8 =  Eq(-Integral(-1/(-(1 - (1 - _u2**2)**(1/2))*_u2 + _u2), (_u2, __a, \
    x/f(x))) + log(C1*f(x)), 0)
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3)
    assert checksol(eq4, f(x), sol4)
    assert checksol(eq5, f(x), sol5)
    assert checksol(eq6, f(x), sol6)
    assert checksol(eq7, f(x), sol7)
    assert checksol(eq8, f(x), sol8)

def test_homogeneous_coeff_ode2():
    eq1 = f(x).diff(x) - f(x)/x+1/sin(f(x)/x)
    eq2 = x**2 + f(x)**2 - 2*x*f(x)*f(x).diff(x)
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol1 = Eq(f(x), x*acos(log(C1*x)))
    sol2 = [Eq(f(x), sqrt(C1*x + x**2)), Eq(f(x), -sqrt(C1*x + x**2))]
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2[0])
    assert checksol(eq2, f(x), sol2[1])
    assert checksol(eq3, f(x), sol3)

@XFAIL
def test_homogeneous_coeff_ode3():
    # This fails because of issue 1429
    eq = x+f(x)-(x-f(x))*f(x).diff(x)
    sol = Eq(f(x), 2*x*atan(C1*x))
    assert dsolve(eq, f(x)) == sol
    assert checksol(eq, f(x), sol)

def test_homogeneous_coeff_ode4():
    # This can be solved explicitly, but the the integration engine cannot handle
    # it (see issue 1452).  The explicit solution is included in an XFAIL test
    # below. checksol fails for this equation, so its test is in
    # test_homogeneous_order_ode1_sol above. It has to compare string
    # expressions because u2 and _a are dummy variables.
    eq = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    solstr =  "-Integral(-1/(-(1 - (1 - _u2**2)**(1/2))*_u2 + _u2), (_u2, " + \
    "__a, x/f(x))) + log(C1*f(x)) == 0"
    assert str(dsolve(eq, f(x))) == solstr


@XFAIL
def test_homogeneous_order_ode4_explicit():
    eq = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    sol = Eq(f(x)**2-C1*x, f(x)*sqrt(f(x)**2-x**2))
    assert dsolve(eq, f(x)) == sol
    assert checksol(eq, f(x)) == sol

def test_homogeneous_norder():
    # From Exercise 20, in Ordinary Differential Equations, Tenenbaum and Pollard
    a = Symbol('a', positive=True)
    k = Symbol('k', real=True)
    eq1 = f(x).diff(x, 2) + 2*f(x).diff(x)
    eq2 = f(x).diff(x, 2) - 3*f(x).diff(x) + 2*f(x)
    eq3 = f(x).diff(x, 2) - f(x)
    eq4 = f(x).diff(x, 3) + f(x).diff(x, 2) - 6*f(x).diff(x)
    eq5 = 6*f(x).diff(x, 2) - 11*f(x).diff(x) + 4*f(x)
    eq6 = Eq(f(x).diff(x, 2) + 2*f(x).diff(x) - f(x), 0)
    eq7 = diff(f(x), x, 3) + diff(f(x), x, 2) - 10*diff(f(x), x) - 6*f(x)
    eq8 = f(x).diff(x, 4) - f(x).diff(x, 3) - 4*f(x).diff(x, 2) + 4*f(x).diff(x)
    eq9 = f(x).diff(x, 4) + 4*f(x).diff(x, 3) + f(x).diff(x, 2) - \
        4*f(x).diff(x) - 2*f(x)
    eq10 = f(x).diff(x, 4) - a**2*f(x)
    eq11 = f(x).diff(x, 2) - 2*k*f(x).diff(x) - 2*f(x)
    eq12 = f(x).diff(x, 2) + 4*k*f(x).diff(x) - 12*k**2*f(x)
    eq13 = f(x).diff(x, 4)
    eq14 = f(x).diff(x, 2) + 4*f(x).diff(x) + 4*f(x)
    eq15 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x)
    eq16 = f(x).diff(x, 3) - 6*f(x).diff(x, 2) +12*f(x).diff(x) - 8*f(x)
    eq17 = f(x).diff(x, 2) - 2*a*f(x).diff(x) + a**2*f(x)
    eq18 = f(x).diff(x, 4) + 3*f(x).diff(x, 3)
    eq19 = f(x).diff(x, 4) - 2*f(x).diff(x, 2)
    eq20 = f(x).diff(x, 4) + 2*f(x).diff(x, 3) - 11*f(x).diff(x, 2) - \
        12*f(x).diff(x) + 36*f(x)
    eq21 = 36*f(x).diff(x, 4) - 37*f(x).diff(x, 2) + 4*f(x).diff(x) + 5*f(x)
    eq22 = f(x).diff(x, 4) - 8*f(x).diff(x, 2) + 16*f(x)
    eq23 = f(x).diff(x, 2) - 2*f(x).diff(x) + 5*f(x)
    eq24 = f(x).diff(x, 2) - f(x).diff(x) + f(x)
    eq25 = f(x).diff(x, 4) + 5*f(x).diff(x, 2) + 6*f(x)
    eq26 = f(x).diff(x, 2) - 4*f(x).diff(x) + 20*f(x)
    eq27 = f(x).diff(x, 4) + 4*f(x).diff(x, 2) + 4*f(x)
    eq28 = f(x).diff(x, 3) + 8*f(x)
    eq29 = f(x).diff(x, 4) + 4*f(x).diff(x, 2)
    eq30 = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x)
    sol1 = Eq(f(x), C1 + C2*exp(-2*x))
    sol2 = Eq(f(x), C1*exp(x) + C2*exp(2*x))
    sol3 = Eq(f(x), C1*exp(x) + C2*exp(-x))
    sol4 = Eq(f(x), C1 + C2*exp(-3*x) + C3*exp(2*x))
    sol5 = Eq(f(x), C1*exp(x/2) + C2*exp(4*x/3))
    sol6 = Eq(f(x), C1*exp(-x*(1 + sqrt(2))) + C2*exp(-x*(1 - sqrt(2))))
    sol7 = Eq(f(x), C1*exp(3*x) + C2*exp(-x*(2 + sqrt(2))) + C3*exp(-x*(2 - sqrt(2))))
    sol8 = Eq(f(x), C1 + C2*exp(x) + C3*exp(-2*x) + C4*exp(2*x))
    sol9 = Eq(f(x), C1*exp(x) + C2*exp(-x) + C3*exp(-x*(2 + sqrt(2))) + C4*exp(-x*(2 - \
        sqrt(2))))
    sol10 = Eq(f(x), C1*sin(x*sqrt(a)) + C2*cos(x*sqrt(a)) + C3*exp(x*sqrt(a)) + \
        C4*exp(-x*sqrt(a)))
    sol11 = Eq(f(x), C1*exp(x*(k + sqrt(8 + 4*k**2)/2)) + C2*exp(x*(k - sqrt(8 + \
        4*k**2)/2)))
    sol12 = Eq(f(x), C1*exp(x*(-4*abs(k) - 2*k)) + C2*exp(x*(-2*k + 4*abs(k))))
    sol13 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*x**3)
    sol14 = Eq(f(x), (C1 + C2*x)*exp(-2*x))
    sol15 = Eq(f(x), (C1 + C2*x)*exp(-x) + C3*exp(x/3))
    sol16 = Eq(f(x), (C1 + C2*x + C3*x**2)*exp(2*x))
    sol17 = Eq(f(x), (C1 + C2*x)*exp(a*x))
    sol18 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*exp(-3*x))
    sol19 = Eq(f(x), C1 + C2*x + C3*exp(x*sqrt(2)) + C4*exp(-x*sqrt(2)))
    sol20 = Eq(f(x), (C1 + C2*x)*exp(-3*x) + (C3 + C4*x)*exp(2*x))
    sol21 = Eq(f(x), C1*exp(x/2) + C2*exp(-x) + C3*exp(-x/3) + C4*exp(5*x/6))
    sol22 = Eq(f(x), (C1 + C2*x)*exp(-2*x) + (C3 + C4*x)*exp(2*x))
    sol23 = Eq(f(x), (C1*sin(2*x) + C2*cos(2*x))*exp(x))
    sol24 = Eq(f(x), (C1*sin(x*sqrt(3)/2) + C2*cos(x*sqrt(3)/2))*exp(x/2))
    sol25 = Eq(f(x), C1*cos(x*sqrt(3)) + C2*sin(x*sqrt(3)) + C3*sin(x*sqrt(2)) + C4*cos(x*sqrt(2)))
    sol26 = Eq(f(x), (C1*sin(4*x) + C2*cos(4*x))*exp(2*x))
    sol27 = Eq(f(x), (C1 + C2*x)*sin(x*sqrt(2)) + (C3 + C4*x)*cos(x*sqrt(2)))
    sol28 = Eq(f(x), (C1*sin(x*sqrt(3)) + C2*cos(x*sqrt(3)))*exp(x) + C3*exp(-2*x))
    sol29 = Eq(f(x), C1 + C2*sin(2*x) + C3*cos(2*x) + C4*x)
    sol30 = Eq(f(x), C1 + (C2 + C3*x)*sin(x) + (C4 + C5*x)*cos(x))
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert dsolve(eq4, f(x)) == sol4
    assert dsolve(eq5, f(x)) == sol5
    assert dsolve(eq6, f(x)) == sol6
    assert dsolve(eq7, f(x)) == sol7
    assert dsolve(eq8, f(x)) == sol8
    assert dsolve(eq9, f(x)) == sol9
    assert dsolve(eq10, f(x)) == sol10
    assert dsolve(eq11, f(x)) == sol11
    assert dsolve(eq12, f(x)) == sol12
    assert dsolve(eq13, f(x)) == sol13
    assert dsolve(eq14, f(x)) == sol14
    assert dsolve(eq15, f(x)) == sol15
    assert dsolve(eq16, f(x)) == sol16
    assert dsolve(eq17, f(x)) == sol17
    assert dsolve(eq18, f(x)) == sol18
    assert dsolve(eq19, f(x)) == sol19
    assert dsolve(eq20, f(x)) == sol20
    assert dsolve(eq21, f(x)) == sol21
    assert dsolve(eq22, f(x)) == sol22
    assert dsolve(eq23, f(x)) == sol23
    assert dsolve(eq24, f(x)) == sol24
    assert dsolve(eq25, f(x)) == sol25
    assert dsolve(eq26, f(x)) == sol26
    assert dsolve(eq27, f(x)) == sol27
    assert dsolve(eq28, f(x)) == sol28
    assert dsolve(eq29, f(x)) == sol29
    assert dsolve(eq30, f(x)) == sol30
    assert checksol(eq1, f(x), sol1)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3)
    assert checksol(eq4, f(x), sol4)
    assert checksol(eq5, f(x), sol5)
    assert checksol(eq6, f(x), sol6)
    assert checksol(eq7, f(x), sol7)
    assert checksol(eq8, f(x), sol8)
    assert checksol(eq9, f(x), sol9)
    assert checksol(eq10, f(x), sol10)
    assert checksol(eq11, f(x), sol11)
    assert checksol(eq12, f(x), sol12)
    assert checksol(eq13, f(x), sol13)
    assert checksol(eq14, f(x), sol14)
    assert checksol(eq15, f(x), sol15)
    assert checksol(eq16, f(x), sol16)
    assert checksol(eq17, f(x), sol17)
    assert checksol(eq18, f(x), sol18)
    assert checksol(eq19, f(x), sol19)
    assert checksol(eq20, f(x), sol20)
    assert checksol(eq21, f(x), sol21)
    assert checksol(eq22, f(x), sol22)
    assert checksol(eq23, f(x), sol23)
    assert checksol(eq24, f(x), sol24)
    assert checksol(eq25, f(x), sol25)
    assert checksol(eq26, f(x), sol26)
    assert checksol(eq27, f(x), sol27)
    assert checksol(eq28, f(x), sol28)
    assert checksol(eq29, f(x), sol29)
    assert checksol(eq30, f(x), sol30)

def test_Liouville_ODE():
    # First part used to be test_ODE_1() from test_solvers.py
    eq1 = diff(f(x),x)/x+diff(f(x),x,x)/2- diff(f(x),x)**2/2
    eq1a = diff(x*exp(-f(x)), x, x)
    eq2 = eq1*exp(-f(x))/exp(f(x))
    eq3 = diff(f(x), x, x) + 1/f(x)*(diff(f(x), x))**2 + 1/x*diff(f(x), x)
    eq4 = x*diff(f(x), x, x) + x/f(x)*diff(f(x), x)**2 + x*diff(f(x), x)
    eq5 = Eq((x*exp(f(x))).diff(x, x), 0)
    sol1 = Eq(C1 + C2/x - exp(-f(x)), 0)
    # If solve() is ever improved, this is a better solution
    sol1a = Eq(f(x), -log((C1*x+C2)/x))
    sol2 = Eq(f(x), -log(C1 + C2/x))
    sol3 = [Eq(f(x), -sqrt(C1 + C2*log(x))), Eq(f(x), sqrt(C1 + C2*log(x)))]
    sol4 = [Eq(f(x), sqrt(C1 + C2*exp(-x))), Eq(f(x), -sqrt(C1 + C2*exp(-x)))]
    sol5 = Eq(f(x), -log(x) + log(C1 + C2*x))
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq1a, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert dsolve(eq4, f(x)) == sol4
    assert dsolve(eq5, f(x)) == sol5
    assert checksol(sol1, f(x), sol1a)
    assert checksol(eq1, f(x), sol1a)
    assert checksol(eq1a, f(x), sol1a)
    assert checksol(eq2, f(x), sol2)
    assert checksol(eq3, f(x), sol3[0])
    assert checksol(eq3, f(x), sol3[1])
    assert checksol(eq4, f(x), sol4[0])
    assert checksol(eq4, f(x), sol4[1])
    assert checksol(eq5, f(x), sol5)

