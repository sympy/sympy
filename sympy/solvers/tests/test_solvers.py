from sympy import Matrix, Symbol, solve, exp, log, cos, acos, Rational, Eq, \
        sqrt, oo, LambertW, pi, I, sin, asin, Function, diff, Derivative, \
        symbols, S, raises, sympify, var, simplify
from sympy.solvers import solve_linear_system, solve_linear_system_LU,dsolve,\
     tsolve, deriv_degree

from sympy.solvers.solvers import guess_solve_strategy, GS_POLY, GS_POLY_CV_1, GS_TRASCENDENTAL, \
        GS_RATIONAL, GS_RATIONAL_CV_1

from sympy.utilities.pytest import XFAIL

def test_guess_strategy():
    """
    See solvers._guess_solve_strategy
    """
    x, y = symbols('xy')

    # polynomial equations
    assert guess_solve_strategy( x, x ) == GS_POLY
    assert guess_solve_strategy( 2*x, x ) == GS_POLY
    assert guess_solve_strategy( x + sqrt(2), x) == GS_POLY
    assert guess_solve_strategy( x + 2**Rational(1,4), x) == GS_POLY
    assert guess_solve_strategy( x**2 + 1, x ) == GS_POLY
    assert guess_solve_strategy( x**2 - 1, x ) == GS_POLY
    assert guess_solve_strategy( x*y + y, x ) == GS_POLY
    assert guess_solve_strategy( x*exp(y) + y, x) == GS_POLY
    assert guess_solve_strategy( (x - y**3)/(y**2*(1 - y**2)**(1/2)), x) == GS_POLY

    # polynomial equations via a change of variable
    assert guess_solve_strategy( x**Rational(1,2) + 1, x ) == GS_POLY_CV_1
    assert guess_solve_strategy( x**Rational(1,3) + x**Rational(1,2) + 1, x ) == GS_POLY_CV_1

    # polynomial equation multiplying both sides by x**n
    assert guess_solve_strategy( x + 1/x + y, x )

    # rational functions
    assert guess_solve_strategy( (x+1)/(x**2 + 2), x) == GS_RATIONAL
    assert guess_solve_strategy( (x - y**3)/(y**2*(1 - y**2)**(1/2)), y) == GS_RATIONAL

    # rational functions via the change of variable y -> x**n
    assert guess_solve_strategy( (x**Rational(1,2) + 1)/(x**Rational(1,3) + x**Rational(1,2) + 1), x ) \
                                == GS_RATIONAL_CV_1

    #trascendental functions
    assert guess_solve_strategy( exp(x) + 1, x ) == GS_TRASCENDENTAL
    assert guess_solve_strategy( 2*cos(x)-y, x ) == GS_TRASCENDENTAL
    assert guess_solve_strategy( exp(x) + exp(-x) - y, x ) == GS_TRASCENDENTAL

def test_solve_polynomial():
    x, y = map(Symbol, 'xy')


    assert solve(3*x-2, x) == [Rational(2,3)]
    assert solve(Eq(3*x, 2), x) == [Rational(2,3)]

    assert solve(x**2-1, x) in [[-1, 1], [1, -1]]
    assert solve(Eq(x**2, 1), x) in [[-1, 1], [1, -1]]

    assert solve( x - y**3, x) == [y**3]
    assert solve( x - y**3, y) in [
    [
            (-x**Rational(1,3))/2 + I*sqrt(3)*x**Rational(1,3)/2,
            x**Rational(1,3),
            (-x**Rational(1,3))/2 - I*sqrt(3)*x**Rational(1,3)/2
    ],
    [
        (-x**Rational(1,3))/2 + I*sqrt(3)*x**Rational(1,3)/2,
        (-x**Rational(1,3))/2 - I*sqrt(3)*x**Rational(1,3)/2,
        x**Rational(1,3)
    ]
    ]


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

def test_solve_polynomial_cv_1():
    """
    Test for solving on equations that can be converted to a polynomial equation
    using the change of variable y -> x**Rational(p, q)
    """

    x = Symbol('x')

    assert solve( x**Rational(1,2) - 1, x) == [1]
    assert solve( x**Rational(1,2) - 2, x) == [sqrt(2)]

def test_solve_polynomial_cv_2():
    """
    Test for solving on equations that can be converted to a polynomial equation
    multiplying both sides of the equation by x**m
    """

    x = Symbol('x')

    assert solve( x + 1/x - 1, x) == [Rational(1,2) + I*sqrt(3)/2, Rational(1,2) - I*sqrt(3)/2]

def test_solve_rational():
    x = Symbol('x')
    y = Symbol('y')

    solve( ( x - y**3 )/( (y**2)*sqrt(1 - y**2) ), x) == [x**Rational(1,3)]

def test_linear_system():
    x, y, z, t, n = map(Symbol, 'xyztn')

    assert solve([x-1, x-y, x-2*y, y-1], [x,y]) is None

    assert solve([x-1, x-y, x-2*y, x-1], [x,y]) is None
    assert solve([x-1, x-1, x-y, x-2*y], [x,y]) is None

    assert solve([x+5*y-2, -3*x+6*y-15], x, y) == {x: -3, y: 1}

    M = Matrix([[0,0,n*(n+1),(n+1)**2,0],
                [n+1,n+1,-2*n-1,-(n+1),0],
                [-1, 0, 1, 0, 0]])

    assert solve_linear_system(M, x, y, z, t) == \
           {y: 0, z: (-t-t*n)/n, x: (-t-t*n)/n}

def test_linear_systemLU():
    x, y, z, n = map(Symbol, 'xyzn')

    M = Matrix([[1,2,0,1],[1,3,2*n,1],[4,-1,n**2,1]])

    assert solve_linear_system_LU(M, [x,y,z]) == {z: -3/(n**2+18*n),
                                                  x: 1-12*n/(n**2+18*n),
                                                  y: 6*n/(n**2+18*n)}

def test_ODE_first_order():
    f = Function('f')
    x = Symbol('x')
    assert dsolve(3*f(x).diff(x) -1, f(x)) == x/3 + Symbol("C1")
    assert dsolve(x*f(x).diff(x) -1, f(x)) == log(x) + Symbol("C1")

def test_ODE_second_order():
    f = Function('f')
    x, C1, C2 = map(Symbol, ['x', 'C1', 'C2'])
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
def test_tsolve():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    assert solve(exp(x)-3, x) == [log(3)]
    assert solve(cos(x)-y, x) == [acos(y)]
    assert solve(2*cos(x)-y,x)== [acos(y/2)]
    # XXX in the following test, log(2*y + 2*...) should -> log(2) + log(y +...)
    assert solve(exp(x)+exp(-x)-y,x)    == [-log(4) + log(2*y + 2*(-4 + y**2)**Rational(1,2)),
                                            -log(4) + log(2*y - 2*(-4 + y**2)**Rational(1,2))]
    assert tsolve(exp(x)-3, x) == [log(3)]
    assert tsolve(Eq(exp(x), 3), x) == [log(3)]
    assert tsolve(log(x)-3, x) == [exp(3)]
    assert tsolve(sqrt(3*x)-4, x) == [Rational(16,3)]
    assert tsolve(3**(x+2), x) == [-oo]
    assert tsolve(3**(2-x), x) == [oo]
    assert tsolve(4*3**(5*x+2)-7, x) == [(log(Rational(7,4))-2*log(3))/(5*log(3))]
    assert tsolve(x+2**x, x) == [-LambertW(log(2))/log(2)]
    assert tsolve(3*x+5+2**(-5*x+3), x) == \
        [-Rational(5,3) + LambertW(-10240*2**Rational(1,3)*log(2)/3)/(5*log(2))]
    assert tsolve(5*x-1+3*exp(2-7*x), x) == \
        [Rational(1,5) + LambertW(-21*exp(Rational(3,5))/5)/7]
    assert tsolve(2*x+5+log(3*x-2), x) == \
        [Rational(2,3) + LambertW(2*exp(-Rational(19,3))/3)/2]
    assert tsolve(3*x+log(4*x), x) == [LambertW(Rational(3,4))/3]
    assert tsolve((2*x+8)*(8+exp(x)), x) == [-4]
    assert tsolve(2*exp(3*x+4)-3, x) == [-Rational(4,3)+log(Rational(3,2))/3]
    assert tsolve(2*log(3*x+4)-3, x) == [(exp(Rational(3,2))-4)/3]
    assert tsolve(exp(x)+1, x) == [pi*I]
    assert tsolve(x**2 - 2**x, x) == [2]
    assert tsolve(x**3 - 3**x, x) == [-3/log(3)*LambertW(-log(3)/3)]
    assert tsolve(2*(3*x+4)**5 - 6*7**(3*x+9), x) == \
        [Rational(-4,3) - 5/log(7)/3*LambertW(-7*2**Rational(4,5)*6**Rational(1,5)*log(7)/10)]

    assert tsolve(z*cos(x)-y, x)      == [acos(y/z)]
    assert tsolve(z*cos(2*x)-y, x)    == [acos(y/z)/2]
    assert tsolve(z*cos(sin(x))-y, x) == [asin(acos(y/z))]

    assert tsolve(z*cos(x), x)        == [acos(0)]

    assert tsolve(exp(x)+exp(-x)-y, x)== [log(y/2 + Rational(1,2)*(-4 + y**2)**Rational(1,2)),
                                          log(y/2 - Rational(1,2)*(-4 + y**2)**Rational(1,2))]
