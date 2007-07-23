
from sympy import *

from sympy.modules.matrices import Matrix
from sympy.modules.solvers import solve_linear_system, solve_linear_system_LU

def test_solve():
    x, y = map(Symbol, 'xy')

    assert solve(3*x-2, x) == [Rational(2,3)]
    assert solve(3*x == 2, x) == [Rational(2,3)]

    assert solve(x**2-1, x) == [-1, 1]
    assert solve(x**2 == 1, x) == [-1, 1]

def _test_linear_system(): # needs simplify()
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
