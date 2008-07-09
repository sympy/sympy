from sympy.solvers.numeric import maxnorm
from sympy.solvers import msolve
from sympy.utilities.lambdify import lambdify
from sympy import Symbol, Matrix

def test_msolve():
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    f1 = 3 * x1**2 - 2 * x2**2 - 1
    f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
    f = Matrix(f1, f2).T
    F = lambdify((x1, x2), f.T)
    # numeric.newton is tested in this example too
    for x0 in ((-1., 1.), (1., -2.), (4., 4.), (-4., -4.)):
        x = msolve((x1, x2), f, x0, tol=1.e-8)
        assert maxnorm(F(*x)) <= 1.e-11
    # TODO: 3-dimensional system
