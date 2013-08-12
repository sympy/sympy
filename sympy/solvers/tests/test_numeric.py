from sympy import Eq, Matrix, pi, sin, sqrt, Symbol, Integral, Piecewise, symbols
from sympy.mpmath import mnorm, mpf
from sympy.solvers import nsolve
from sympy.utilities.lambdify import lambdify
from sympy.utilities.pytest import raises, XFAIL


def test_nsolve():
    # onedimensional
    x = Symbol('x')
    assert nsolve(sin(x), 2) - pi.evalf() < 1e-15
    assert nsolve(Eq(2*x, 2), x, -10) == nsolve(2*x - 2, -10)
    # Testing checks on number of inputs
    raises(TypeError, lambda: nsolve(Eq(2*x, 2)))
    raises(TypeError, lambda: nsolve(Eq(2*x, 2), x, 1, 2))
    # Issue 1730
    assert nsolve(x**2/(1 - x)/(1 - 2*x)**2 - 100, x, 0)  # doesn't fail
    # multidimensional
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    f1 = 3 * x1**2 - 2 * x2**2 - 1
    f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
    f = Matrix((f1, f2)).T
    F = lambdify((x1, x2), f.T, modules='mpmath')
    for x0 in [(-1, 1), (1, -2), (4, 4), (-4, -4)]:
        x = nsolve(f, (x1, x2), x0, tol=1.e-8)
        assert mnorm(F(*x), 1) <= 1.e-10
    # The Chinese mathematician Zhu Shijie was the very first to solve this
    # nonlinear system 700 years ago (z was added to make it 3-dimensional)
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    f1 = -x + 2*y
    f2 = (x**2 + x*(y**2 - 2) - 4*y) / (x + 4)
    f3 = sqrt(x**2 + y**2)*z
    f = Matrix((f1, f2, f3)).T
    F = lambdify((x, y, z), f.T, modules='mpmath')

    def getroot(x0):
        root = nsolve(f, (x, y, z), x0)
        assert mnorm(F(*root), 1) <= 1.e-8
        return root
    assert list(map(round, getroot((1, 1, 1)))) == [2.0, 1.0, 0.0]
    assert nsolve([Eq(
        f1), Eq(f2), Eq(f3)], [x, y, z], (1, 1, 1))  # just see that it works
    a = Symbol('a')
    assert nsolve(1/(0.001 + a)**3 - 6/(0.9 - a)**3, a, 0.3).ae(
        mpf('0.31883011387318591'))


def test_issue_3309():
    x = Symbol('x')
    assert nsolve(Piecewise((x, x < 1), (x**2, True)), x, 2) == 0.0


@XFAIL
def test_issue_3309_fail():
    x, y = symbols('x y')
    assert nsolve(Integral(x*y, (x, 0, 5)), y, 2) == 0.0
