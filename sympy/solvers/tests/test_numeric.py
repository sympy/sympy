from sympy import Eq, Expr, Matrix, pi, sin, sqrt, Symbol, tanh, Tuple
from sympy.mpmath import mnorm, mpf
from sympy.solvers import nsolve
from sympy.utilities.lambdify import lambdify
from sympy.utilities.pytest import raises

def test_nsolve():
    # onedimensional
    x = Symbol('x')
    assert nsolve(sin(x), 2) - pi.evalf() < 1e-15
    assert nsolve(Eq(2*x, 2), x, -10) == nsolve(2*x - 2, -10)
    # Testing checks on number of inputs
    raises(TypeError, "nsolve(Eq(2*x,2))")
    raises(TypeError, "nsolve(Eq(2*x,2),x,1,2)")
    # Issue 1730
    assert nsolve(x**2/(1-x)/(1-2*x)**2-100, x, 0) # doesn't fail
    # multidimensional
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    f1 = 3 * x1**2 - 2 * x2**2 - 1
    f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
    f = Matrix((f1, f2)).T
    F = lambdify((x1, x2), f.T, modules='mpmath')
    for x0 in [(-1, 1), (1, -2), (4, 4), (-4, -4)]:
        x = nsolve(f, (x1, x2), x0, tol=1.e-8)
        assert mnorm(F(*x),1) <= 1.e-10
    # The Chinese mathematician Zhu Shijie was the very first to solve this
    # nonlinear system 700 years ago (z was added to make it 3-dimensional)
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    f1 = -x + 2*y
    f2 = (x**2 + x*(y**2 - 2) - 4*y)  /  (x + 4)
    f3 = sqrt(x**2 + y**2)*z
    f = Matrix((f1, f2, f3)).T
    F = lambdify((x, y, z), f.T, modules='mpmath')
    def getroot(x0):
        root = nsolve(f, (x, y, z), x0)
        assert mnorm(F(*root),1) <= 1.e-8
        return root
    assert map(round, getroot((1, 1, 1))) == [2.0, 1.0, 0.0]
    assert nsolve([Eq(f1), Eq(f2), Eq(f3)], [x, y, z], (1, 1, 1)) # just see that it works
    a = Symbol('a')
    assert nsolve(1/(0.001 + a)**3 - 6/(0.9 - a)**3, a, 0.3).ae(
        mpf('0.31883011387318591'))

def test_nsolve_symbolic_expressions():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    raises(TypeError, 'nsolve(tanh(x*y)-x,10)')

    r1 = nsolve(tanh(x*y)-x,x,10)
    assert isinstance(r1, Expr)
    assert r1.evalf(subs={y : 0.2}) == 0

    r2 = nsolve((y-x, x+z), (x, z), (0, 0), imp_function='root2')
    assert isinstance(r2, Expr)
    # TODO: evalf does not play well with mpmath matrix.
    # For the moment I'm only testing _imp_.
    #assert r2.evalf(subs={y : 1.2})[0] == 1.2
    assert r2._imp_(1.2)[0] == 1.2
    #assert r2.evalf(subs={y : 1.2})[1] == -1.2
    assert r2._imp_(1.2)[1] == -1.2

    # testing all combinations of Tuple,list,Matrix arguments
    # TODO: fix evalf so you can use it directly and not _imp_
    # TODO: fix the commented out tests
    t = Tuple(x-y,z-y)
    l = list(t)
    m = Matrix(l)
    assert nsolve(t,(x,z),(0,0))._imp_(1)[0] == 1.0
    assert nsolve(l,(x,z),(0,0))._imp_(1)[0] == 1.0
    assert nsolve(m,(x,z),(0,0))._imp_(1)[0] == 1.0
    assert nsolve(t,(x,z),(0,0))._imp_(1)[0] == 1.0
    assert nsolve(t,[x,z],(0,0))._imp_(1)[0] == 1.0
    #assert nsolve(t,Matrix([x,z]),(0,0))._imp_(1)[0] == 1.0
    # it return nsolve_root(z,y)
    assert nsolve(t,(x,z),[0,0])._imp_(1)[0] == 1.0
    #assert nsolve(t,(x,z),Matrix([0,0]))._imp_(1)[0] == 1.0
    # some strange error: TypeError: cannot create mpf from [0]
    # AND returns [0]

    # testing the arg_order argument
    assert (nsolve(x*y/z-1.,x,10.,arg_order=[y,z])._imp_(10.,1. )
         == nsolve(x*y/z-1.,x,10.,arg_order=[z,y])._imp_(1. ,10.))
