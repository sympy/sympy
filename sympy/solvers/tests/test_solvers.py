
from sympy import *
from sympy.utilities.pytest import XFAIL

from sympy.matrices import Matrix
from sympy.solvers import solve_linear_system, solve_linear_system_LU,dsolve

def test_solve():
    x, y = map(Symbol, 'xy')

    assert solve(3*x-2, x) == [Rational(2,3)]
    assert solve(3*x == 2, x) == [Rational(2,3)]

    assert solve(x**2-1, x) in [[-1, 1], [1, -1]]
    assert solve(x**2 == 1, x) in [[-1, 1], [1, -1]]

@XFAIL
def test_linear_system(): # needs simplify()
    x, y, z, t, n = map(Symbol, 'xyztn')

    assert solve([x+5*y-2, -3*x+6*y-15], [x, y]) == {x: -3, y: 1}

    M = Matrix( [0,0,n*(n+1),(n+1)**2,0],
                [n+1,n+1,-2*n-1,-(n+1),0],
                [-1, 0, 1, 0, 0] )

    assert solve_linear_system(M, [x, y, z, t]) == \
           {y: 0, z: (-t-t*n)/n, x: (-t-t*n)/n}

def test_linear_systemLU():
    x, y, z, n = map(Symbol, 'xyzn')

    M = Matrix( [1,2,0,1],[1,3,2*n,1],[4,-1,n**2,1])

    assert solve_linear_system_LU(M, [x,y,z]) == {z: -3/(n**2+18*n),
                                                  x: 1-12*n/(n**2+18*n),
                                                  y: 6*n/(n**2+18*n)}

@XFAIL
def test_ODE_first_order():
    f = Function2('f')
    x = Symbol('x')
    assert dsolve(3*f(x).diff(x) -1, f(x)) == x/3 + Symbol("C1")
    assert dsolve(x*f(x).diff(x) -1, f(x)) == log(x) + Symbol("C1")

@XFAIL
def test_ODE_second_order():
    f = Function2('f')
    x, C1, C2 = map(Symbol, ['x', 'C1', 'C2'])
    assert dsolve(Derivative(f(x),x,x) + 9*f(x), [f(x)]) in \
        [sin(3*x)*C1 + cos(3*x)*C2, sin(3*x)*C2 + cos(3*x)*C1]

@XFAIL
def test_ODE_1():
    l = Function2('l')
    r = Symbol('r')

    e = Derivative(l(r),r)/r+Derivative(l(r),r,r)/2- \
        Derivative(l(r),r)**2/2
    sol = dsolve(e, [l(r)])
    assert (e.subs(l(r), sol)).expand() == 0

    e = e*exp(-l(r))/exp(l(r))
    sol = dsolve(e, [l(r)])
    assert (e.subs(l(r), sol)).expand() == 0
