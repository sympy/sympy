from sympy import Function, dsolve, Symbol, sin, cos, sinh, acos, tan, cosh, \
        I, exp, log, simplify, normal, together, ratsimp, powsimp, \
        fraction, radsimp, Eq, sqrt, pi, erf,  diff, Rational, asinh, trigsimp, \
        S, RootOf, Poly, Integral, atan, Equality, solve, O, LambertW
from sympy.abc import x, y, z
from sympy.solvers.ide import classify_ide 
from sympy.utilities.pytest import XFAIL, skip, raises

n= Symbol('n')
a = Symbol('a')
b = Symbol('b')
C1 = Symbol('C1')
C2 = Symbol('C2')

f = Function('f')
g = Function('g')
K = function('K')

vfnh = ['Non-homogenous', 'Volterra', 'First Kind']
vsnh = ['Non-homogenous', 'Volterra', 'Second Kind']
ffnh = ['Non-homogenous', 'Fredholm', 'First Kind']
fsnh = ['Non-homogenous', 'Fredholm', 'Second Kind']

def test_checkodesol():
    # For the most part, checkodesol is well tested in the tests below.
    # These tests only handle cases not checked below.
    raises(ValueError, "checkodesol(f(x).diff(x), f(x), x)")
    raises(ValueError, "checkodesol(f(x).diff(x), f(x, y), Eq(f(x), x))")
    
def test_classify_ode():
    """
    Tests for volterra integral equation of the first kind
    """
    eq1 = Eq(Integral(f(y)*(x-y),(y,0,x)),g(x))
    assert classify_ide(eq1,f(x)) == vfnh
    eq2 = Eq(Integral(f(y)*(x-y)**n,(y,0,x)),g(x))
    assert classify_ide(eq2,f(x)) == vfnh
    eq3 = Eq(Integral(f(y)*(exp(x-y)**n+1),(y,0,x)),g(x))
    assert classify_ide(eq3,f(x)) == vfnh
    eq4 = Eq(Integral(f(y)*(log(x-y)**n+1)*(x-y),(y,0,x)),g(x))
    assert classify_ide(eq4,f(x)) == vfnh
    eq5 = Eq(Integral((f(y)*sin(x-y))*cos(x-y),(y,0,x)),g(x))
    assert classify_ide(eq5,f(x)) == vfnh
    """
    Tests for volterra integral equation of the second kind
    """
    eq1 = Eq(f(x) + Integral(f(y)*(x-y),(y,0,x)),g(x))
    assert classify_ide(eq1,f(x)) == vsnh
    eq2 = Eq(f(x) + C2*Integral(f(y)*((x-y)**(-C1)),(y,0,x)),g(x)) # Generalized Abel equation
    assert classify_ide(eq2,f(x)) == vsnh
    eq3 = Eq(n*Integral((f(y)*exp(x-y)*(x-y)),(y,0,x)),g(x) + f(x))
    assert classify_ide(eq3,f(x)) == vsnh
    eq4 = Eq(f(x) + Integral(f(y)*(log(x-y)**n+1)*(x-y),(y,0,x)),g(x))
    assert classify_ide(eq4,f(x)) == vsnh
    eq5 = Eq(Integral((f(y)*K(x-y)),(y,0,x)),g(x) + f(x)) # Renewal equation
    assert classify_ide(eq5,f(x)) == vsnh
    """
    Tests for fredholm integral equation of the first kind
    """
    eq1 = Eq(Integral(f(y)*log(abs(x-y)),(y,a,b)),g(x)) # Carleman's equation
    assert classify_ide(eq1,f(x)) == ffnh
    eq2 = Eq(Integral(f(y)*(x-y)**n,(y,0,b)),g(x))
    assert classify_ide(eq2,f(x)) == ffnh
    eq3 = Eq(Integral(f(y*sin(y))*log(abs(x-y)),(y,0,pi/2)),g(x)) # Schlomilch equation
    assert classify_ide(eq3,f(x)) == ffnh
    eq4 = Eq(Integral(f(y)*cos((x/2-y/2)),(y,0,pi/2)),g(x))
    assert classify_ide(eq4,f(x)) == ffnh
    eq5 = Eq(Integral(f(y)*K(x-y),(y,0,oo)),g(x)) # Wiener-Hopf equation of the first kind
    assert classify_ide(eq5,f(x)) == ffnh
    """
    Tests for fredholm integral equation of the second kind
    """
    eq1 = Eq(f(x) + Integral(f(y)*(x-y),(y,0,1)),g(x))
    assert classify_ide(eq1,f(x)) == fsnh
    eq2 = Eq(f(x) - n*Integral(f(y)*(1/(x-y) + 1/(x+y-2*x*y)),(y,0,1)),g(x)) # Tricomi equation
    assert classify_ide(eq2,f(x)) == fsnh
    eq3 = Eq(n*Integral((f(y)*exp(x-y)*(x-y)),(y,0,)),g(x) + f(x))
    assert classify_ide(eq3,f(x)) == fsnh
    eq4 = Eq(f(x) - Integral(f(y)*K(x-y),(y,0,oo)),g(x)) # Wiener-Hopf equation of the second kind
    assert classify_ide(eq4,f(x)) == fsnh
    eq5 = Eq(f(x) - n*Integral(f(y)*exp(-abs(x-y)),(y,-oo,oo)),g(x)) # Lalesco-Picard equation
    assert classify_ide(eq5,f(x)) == fsnh
    
