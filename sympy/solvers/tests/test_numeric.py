from sympy.solvers.numeric import maxnorm
from sympy.solvers import msolve
from sympy.utilities.lambdify import lambdify
from sympy import Symbol, Matrix,  sqrt

def test_msolve():
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    f1 = 3 * x1**2 - 2 * x2**2 - 1
    f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
    f = Matrix((f1, f2)).T
    F = lambdify((x1, x2), f.T)
    # numeric.newton is tested in this example too
    for x0 in [(-1., 1.), (1., -2.), (4., 4.), (-4., -4.)]:
        x = msolve((x1, x2), f, x0, tol=1.e-8)
        assert maxnorm(F(*x)) <= 1.e-11
    # The Chinese mathematician Zhu Shijie was the very first to solve this
    # nonlinear system 700 years ago (z was added to make it 3-dimensional)
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    f1 = -x + 2*y
    f2 = (x**2 + x*(y**2 - 2) - 4*y)  /  (x + 4)
    f3 = sqrt(x**2 + y**2)*z
    f = Matrix((f1, f2, f3)).T
    F = lambdify((x,  y,  z), f.T)
    def getroot(x0):
        root = msolve((x,  y,  z),  (f1,  f2,  f3),  x0)
        assert maxnorm(F(*root)) <= 1.e-8
        return root
    assert map(round,  getroot((1.,  1.,  1.))) == [2.0,  1.0,  0.0]
