from sympy.mpmath import *
from sympy.mpmath.optimization import *

def test_findroot():
    # old tests, assuming secant
    mp.dps = 15
    assert findroot(lambda x: 4*x-3, mpf(5)).ae(0.75)
    assert findroot(sin, mpf(3)).ae(pi)
    assert findroot(sin, (mpf(3), mpf(3.14))).ae(pi)
    assert findroot(lambda x: x*x+1, mpc(2+2j)).ae(1j)
    # test all solvers with 1 starting point
    f = lambda x: cos(x)
    for solver in [Secant, MNewton, Muller, ANewton]:
        x = findroot(f, 2., solver=solver)
        assert abs(f(x)) < eps
    # test all solvers with interval of 2 points
    for solver in [Secant, Muller, Bisection, Illinois, Pegasus, Anderson,
                   Ridder]:
        x = findroot(f, (1., 2.), solver=solver)
        assert abs(f(x)) < eps
    # test types
    f = lambda x: (x - 2)**2
    assert isinstance(findroot(f, 1, force_type=mpf, tol=1e-10), mpf)
    assert isinstance(findroot(f, 1., force_type=None, tol=1e-10), float)
    assert isinstance(findroot(f, 1, force_type=complex, tol=1e-10), complex)

def test_mnewton():
    f = lambda x: polyval([1,3,3,1],x)
    x = findroot(f, -0.9, solver='mnewton')
    assert abs(f(x)) < eps

def test_anewton():
    f = lambda x: (x - 2)**100
    x = findroot(f, 1., solver=ANewton)
    assert abs(f(x)) < eps

def test_muller():
    f = lambda x: (2 + x)**3 + 2
    x = findroot(f, 1., solver=Muller)
    assert abs(f(x)) < eps

def test_multiplicity():
    for i in xrange(1, 5):
        assert multiplicity(lambda x: (x - 1)**i, 1) == i
    assert multiplicity(lambda x: x**2, 1) == 0
