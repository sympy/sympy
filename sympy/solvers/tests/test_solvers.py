from sympy import Matrix, Symbol, solve, exp, log, cos, acos, Rational, Eq, \
        sqrt, oo, LambertW, pi, I, sin, asin, Function, diff, Derivative, \
        symbols, S, raises, sympify, var, simplify, Integral

from sympy.solvers import solve_linear_system, solve_linear_system_LU,dsolve,\
     tsolve, deriv_degree

from sympy.solvers.solvers import guess_solve_strategy, GS_POLY, GS_POLY_CV_1, GS_POLY_CV_2,\
    GS_TRANSCENDENTAL, GS_RATIONAL, GS_RATIONAL_CV_1

from sympy.utilities.pytest import XFAIL

def test_swap_back():
    x=var('x');f=Function('f',dummy=True)
    assert solve(Eq(log(f(x)), Integral(x, (x, 1, f(x)))), f(x)) == \
    [exp(Integral(x, (x, 1, f(x))))]

def test_guess_poly():
    """
    See solvers.guess_solve_strategy
    """
    x, y, a = symbols('xya')

    # polynomial equations
    assert guess_solve_strategy( S(4), x ) == GS_POLY
    assert guess_solve_strategy( x, x ) == GS_POLY
    assert guess_solve_strategy( x + a, x ) == GS_POLY
    assert guess_solve_strategy( 2*x, x ) == GS_POLY
    assert guess_solve_strategy( x + sqrt(2), x) == GS_POLY
    assert guess_solve_strategy( x + 2**Rational(1,4), x) == GS_POLY
    assert guess_solve_strategy( x**2 + 1, x ) == GS_POLY
    assert guess_solve_strategy( x**2 - 1, x ) == GS_POLY
    assert guess_solve_strategy( x*y + y, x ) == GS_POLY
    assert guess_solve_strategy( x*exp(y) + y, x) == GS_POLY
    assert guess_solve_strategy( (x - y**3)/(y**2*(1 - y**2)**(S(1)/2)), x) == GS_POLY

def test_guess_poly_cv():
    x, y = symbols('xy')
    # polynomial equations via a change of variable
    assert guess_solve_strategy( x**Rational(1,2) + 1, x ) == GS_POLY_CV_1
    assert guess_solve_strategy( x**Rational(1,3) + x**Rational(1,2) + 1, x ) == GS_POLY_CV_1
    assert guess_solve_strategy( 4*x*(1 - sqrt(x)), x ) == GS_POLY_CV_1

    # polynomial equation multiplying both sides by x**n
    assert guess_solve_strategy( x + 1/x + y, x ) == GS_POLY_CV_2

def test_guess_rational_cv():
    # rational functions
    x, y = symbols('xy')
    assert guess_solve_strategy( (x+1)/(x**2 + 2), x) == GS_RATIONAL
    assert guess_solve_strategy( (x - y**3)/(y**2*(1 - y**2)**(S(1)/2)), y) == GS_RATIONAL_CV_1

    # rational functions via the change of variable y -> x**n
    assert guess_solve_strategy( (x**Rational(1,2) + 1)/(x**Rational(1,3) + x**Rational(1,2) + 1), x ) \
                                == GS_RATIONAL_CV_1

def test_guess_transcendental():
    x, y, a, b = symbols('xyab')
    #transcendental functions
    assert guess_solve_strategy( exp(x) + 1, x ) == GS_TRANSCENDENTAL
    assert guess_solve_strategy( 2*cos(x)-y, x ) == GS_TRANSCENDENTAL
    assert guess_solve_strategy( exp(x) + exp(-x) - y, x ) == GS_TRANSCENDENTAL
    assert guess_solve_strategy(3**x-10, x) == GS_TRANSCENDENTAL
    assert guess_solve_strategy(-3**x+10, x) == GS_TRANSCENDENTAL

    assert guess_solve_strategy(a*x**b-y, x) == GS_TRANSCENDENTAL

def test_solve_polynomial1():
    x, y = symbols('xy')

    assert solve(3*x-2, x) == [Rational(2,3)]
    assert solve(Eq(3*x, 2), x) == [Rational(2,3)]

    assert solve(x**2-1, x) in [[-1, 1], [1, -1]]
    assert solve(Eq(x**2, 1), x) in [[-1, 1], [1, -1]]

    assert solve( x - y**3, x) == [y**3]
    assert sorted(solve( x - y**3, y)) == sorted([
        (-x**Rational(1,3))/2 + I*sqrt(3)*x**Rational(1,3)/2,
        x**Rational(1,3),
        (-x**Rational(1,3))/2 - I*sqrt(3)*x**Rational(1,3)/2,
    ])

    a11,a12,a21,a22,b1,b2 = symbols('a11','a12','a21','a22','b1','b2')

    assert solve([a11*x + a12*y - b1, a21*x + a22*y - b2], x, y) == \
        { y : (a11*b2 - a21*b1)/(a11*a22 - a12*a21),
          x : (a22*b1 - a12*b2)/(a11*a22 - a12*a21) }

    solution = {y: S.Zero, x: S.Zero}

    assert solve((x-y, x+y),  x, y ) == solution
    assert solve((x-y, x+y), (x, y)) == solution
    assert solve((x-y, x+y), [x, y]) == solution

    assert solve( x**3 - 15*x - 4, x) == [-2 + 3**Rational(1,2),
                                           4,
                                           -2 - 3**Rational(1,2) ]

    raises(TypeError, "solve(x**2-pi, pi)")
    raises(ValueError, "solve(x**2-pi)")

def test_solve_polynomial2():
    x = Symbol('x')
    assert solve(4, x) == []

def test_solve_polynomial_cv_1a():
    """
    Test for solving on equations that can be converted to a polynomial equation
    using the change of variable y -> x**Rational(p, q)
    """

    x = Symbol('x')
    assert solve( x**Rational(1,2) - 1, x) == [1]
    assert solve( x**Rational(1,2) - 2, x) == [4]
    assert solve( x**Rational(1,4) - 2, x) == [16]
    assert solve( x**Rational(1,3) - 3, x) == [27]

def test_solve_polynomial_cv_1b():
    x, a = symbols('x a')


    assert set(solve(4*x*(1 - a*x**(S(1)/2)), x)) == set([S(0), 1/a**2])
    assert set(solve(x * (x**(S(1)/3) - 3), x)) == set([S(0), S(27)])

def test_solve_polynomial_cv_2():
    """
    Test for solving on equations that can be converted to a polynomial equation
    multiplying both sides of the equation by x**m
    """

    x = Symbol('x')

    assert solve(x + 1/x - 1, x) in \
        [[ Rational(1,2) + I*sqrt(3)/2, Rational(1,2) - I*sqrt(3)/2],
         [ Rational(1,2) - I*sqrt(3)/2, Rational(1,2) + I*sqrt(3)/2]]

def test_solve_rational():
    """Test solve for rational functions"""
    x, y, a, b = symbols('xyab')
    assert solve( ( x - y**3 )/( (y**2)*sqrt(1 - y**2) ), x) == [y**3]
    assert solve(y-b/(1+a*x), x) == [(b - y)/(a*y)]

def test_linear_system():
    x, y, z, t, n = symbols('xyztn')

    assert solve([x-1, x-y, x-2*y, y-1], [x,y]) is None

    assert solve([x-1, x-y, x-2*y, x-1], [x,y]) is None
    assert solve([x-1, x-1, x-y, x-2*y], [x,y]) is None

    assert solve([x+5*y-2, -3*x+6*y-15], x, y) == {x: -3, y: 1}

    M = Matrix([[0,0,n*(n+1),(n+1)**2,0],
                [n+1,n+1,-2*n-1,-(n+1),0],
                [-1, 0, 1, 0, 0]])

    assert solve_linear_system(M, x, y, z, t) == \
           {y: 0, z: -((t+t*n)/n), x: -((t+t*n)/n)}

def test_linear_systemLU():
    x, y, z, n = symbols('xyzn')

    M = Matrix([[1,2,0,1],[1,3,2*n,1],[4,-1,n**2,1]])

    assert solve_linear_system_LU(M, [x,y,z]) == {z: -3/(n**2+18*n),
                                                  x: 1-12*n/(n**2+18*n),
                                                  y: 6*n/(n**2+18*n)}

def test_ODE_first_order():
    f = Function('f')
    x = Symbol('x')
    C1 = Symbol('C1')
    assert dsolve(3*f(x).diff(x) -1, f(x)) == x/3 + C1
    assert dsolve(x*f(x).diff(x) -1, f(x)) == log(x) + C1
    assert dsolve(x*f(x).diff(x)+f(x)-f(x)**2,f(x)) == 1/(x*(C1 + 1/x))
    assert dsolve(cos(f(x))-(x*sin(f(x))-f(x)**2)*f(x).diff(x),f(x)) == \
    Equality(x*cos(f(x))+f(x)**3/3,C1)
    assert dsolve(sin(x)*cos(f(x))+cos(x)*sin(f(x))*f(x).diff(x),f(x)) == \
    Equality(f(x),acos((1-C1)/cos(x)))

def test_ODE_second_order():
    f = Function('f')
    x, C1, C2 = symbols('x C1 C2')
    assert dsolve(Derivative(f(x),x,x) + 9*f(x), [f(x)]) in \
        [sin(3*x)*C1 + cos(3*x)*C2, sin(3*x)*C2 + cos(3*x)*C1]

def test_ODE_1():
    l = Function('l')
    r = Symbol('r')

    e = Derivative(l(r),r)/r+Derivative(l(r),r,r)/2- \
        Derivative(l(r),r)**2/2
    sol = dsolve(e, [l(r)])
    assert (e.subs(l(r), sol)).expand() == 0

    e = e*exp(-l(r))/exp(l(r))
    sol = dsolve(e, [l(r)])
    assert (e.subs(l(r), sol)).expand() == 0

def test_deriv_degree():
    f = Function('f')
    x = Symbol('x')
    assert deriv_degree(3*x*exp(f(x)), f(x)) == 0
    assert deriv_degree(x*diff(f(x),x)+3*x*f(x)-sin(x)/x, f(x)) == 1
    assert deriv_degree(x**2*f(x).diff(x,x)+x*diff(f(x),x)-f(x),f(x)) == 2
    assert deriv_degree(diff(x*exp(f(x)),x,x), f(x)) == 2
    assert deriv_degree(diff(x*diff(x*exp(f(x)), x,x), x), f(x)) == 3

# Note: multiple solutions exist for some of these equations, so the tests
# should be expected to break if the implementation of the solver changes
# in such a way that a different branch is chosen
def test_tsolve_1():
    a, b = symbols('ab')
    x, y, z = symbols('xyz')
    assert solve(exp(x)-3, x) == [log(3)]
    assert solve((a*x+b)*(exp(x)-3), x) == [-b/a, log(3)]
    assert solve(cos(x)-y, x) == [acos(y)]
    assert solve(2*cos(x)-y,x)== [acos(y/2)]
    raises(NotImplementedError, "solve(Eq(cos(x), sin(x)), x)")

    # XXX in the following test, log(2*y + 2*...) should -> log(2) + log(y +...)
    assert solve(exp(x)+exp(-x)-y,x)    == [-log(4) + log(2*y + 2*(-4 + y**2)**Rational(1,2)),
                                            -log(4) + log(2*y - 2*(-4 + y**2)**Rational(1,2))]
    assert solve(exp(x)-3, x) == [log(3)]
    assert solve(Eq(exp(x), 3), x) == [log(3)]
    assert solve(log(x)-3, x) == [exp(3)]
    assert solve(sqrt(3*x)-4, x) == [Rational(16,3)]
    assert solve(3**(x+2), x) == [-oo]
    assert solve(3**(2-x), x) == [oo]
    assert solve(4*3**(5*x+2)-7, x) == [(-log(4) - 2*log(3) + log(7))/(5*log(3))]
    assert solve(x+2**x, x) == [-LambertW(log(2))/log(2)]
    assert solve(3*x+5+2**(-5*x+3), x) in \
        [[-Rational(5,3) + LambertW(-10240*2**Rational(1,3)*log(2)/3)/(5*log(2))],\
        [(-25*log(2) + 3*LambertW(-10240*2**(Rational(1, 3))*log(2)/3))/(15*log(2))]]
    assert solve(5*x-1+3*exp(2-7*x), x) == \
        [Rational(1,5) + LambertW(-21*exp(Rational(3,5))/5)/7]
    assert solve(2*x+5+log(3*x-2), x) == \
        [Rational(2,3) + LambertW(2*exp(-Rational(19,3))/3)/2]
    assert solve(3*x+log(4*x), x) == [LambertW(Rational(3,4))/3]
    assert solve((2*x+8)*(8+exp(x)), x) == [-4, log(8) + pi*I]
    assert solve(2*exp(3*x+4)-3, x) in [ [-Rational(4,3)+log(Rational(3,2))/3],\
                                         [Rational(-4, 3) - log(2)/3 + log(3)/3]]
    assert solve(2*log(3*x+4)-3, x) == [(exp(Rational(3,2))-4)/3]
    assert solve(exp(x)+1, x) == [pi*I]
    assert solve(x**2 - 2**x, x) == [2]
    assert solve(x**3 - 3**x, x) == [-3/log(3)*LambertW(-log(3)/3)]
    assert solve(2*(3*x+4)**5 - 6*7**(3*x+9), x) in \
        [[Rational(-4,3) - 5/log(7)/3*LambertW(-7*2**Rational(4,5)*6**Rational(1,5)*log(7)/10)],\
         [(-5*LambertW(-7*2**(Rational(4, 5))*6**(Rational(1, 5))*log(7)/10) - 4*log(7))/(3*log(7))], \
         [-((4*log(7) + 5*LambertW(-7*2**Rational(4,5)*6**Rational(1,5)*log(7)/10))/(3*log(7)))]]

    assert solve(z*cos(x)-y, x)      == [acos(y/z)]
    assert solve(z*cos(2*x)-y, x)    == [acos(y/z)/2]
    assert solve(z*cos(sin(x))-y, x) == [asin(acos(y/z))]

    assert solve(z*cos(x), x)        == [acos(0)]

    assert solve(exp(x)+exp(-x)-y, x)== [-log(4) + log(2*y + 2*(-4 + y**2)**(Rational(1, 2))),
                                          -log(4) + log(2*y - 2*(-4 + y**2)**(Rational(1, 2)))]
    # issue #1409
    assert solve(y - b*x/(a+x), x) == [a*y/(b - y)]
    assert solve(y - b*exp(a/x), x) == [a/(-log(b) + log(y))]
    # issue #1408
    assert solve(y-b/(1+a*x),x) == [(b - y)/(a*y)]
    # issue #1407
    assert solve(y-a*x**b , x) == [y**(1/b)*(1/a)**(1/b)]
    # issue #1406
    assert solve(z**x - y, x) == [log(y)/log(z)]
    # issue #1405
    assert solve(2**x - 10, x) == [log(10)/log(2)]

def test_tsolve_2():
    x, y, a, b = symbols('xyab')
    assert solve(y-a*x**b, x) == [y**(1/b)*(1/a)**(1/b)]

def test_solve_for_functions_derivatives():
    t = Symbol('t')
    x = Function('x')(t)
    y = Function('y')(t)
    a11,a12,a21,a22,b1,b2 = symbols('a11','a12','a21','a22','b1','b2')

    soln = solve([a11*x + a12*y - b1, a21*x + a22*y - b2], x, y)
    assert soln == { y : (a11*b2 - a21*b1)/(a11*a22 - a12*a21),
        x : (a22*b1 - a12*b2)/(a11*a22 - a12*a21) }

    assert solve(x-1, x) == [1]
    assert solve(3*x-2, x) == [Rational(2,3)]

    soln = solve([a11*x.diff(t) + a12*y.diff(t) - b1, a21*x.diff(t) +
            a22*y.diff(t) - b2], x.diff(t), y.diff(t))
    assert soln == { y.diff(t) : (a11*b2 - a21*b1)/(a11*a22 - a12*a21),
            x.diff(t) : (a22*b1 - a12*b2)/(a11*a22 - a12*a21) }

    assert solve(x.diff(t)-1, x.diff(t)) == [1]
    assert solve(3*x.diff(t)-2, x.diff(t)) == [Rational(2,3)]

    eqns = set((3*x - 1, 2*y-4))
    assert solve(eqns, set((x,y))) == { x : Rational(1, 3), y: 2 }
    x = Symbol('x')
    f = Function('f')
    F = x**2 + f(x)**2 - 4*x - 1
    assert solve(F.diff(x), diff(f(x), x)) == [(2 - x)/f(x)]

    # Mixed cased with a Symbol and a Function
    x = Symbol('x')
    y = Function('y')(t)

    soln = solve([a11*x + a12*y.diff(t) - b1, a21*x +
            a22*y.diff(t) - b2], x, y.diff(t))
    assert soln == { y.diff(t) : (a11*b2 - a21*b1)/(a11*a22 - a12*a21),
            x : (a22*b1 - a12*b2)/(a11*a22 - a12*a21) }


def test_issue626():
    x = Symbol("x")
    f = Function("f")
    F = x**2 + f(x)**2 - 4*x - 1
    e = F.diff(x)
    assert solve(e, f(x).diff(x)) == [(2-x)/f(x)]
