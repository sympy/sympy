from sympy import Function, dsolve, Symbol, sin, cos, sinh, acos, tan, cosh, \
        I, exp, log, simplify, normal, together, ratsimp, powsimp, \
        fraction, radsimp, Eq, sqrt, pi, erf,  diff, Rational, asinh, trigsimp, \
        S, RootOf, Poly, Integral, atan, Equality, solve, O, LambertW
from sympy.abc import x, y, z
from sympy.solvers.ide import classify_ide , solve_series, checkidesol, \
        solve_adomian, solve_approximate 
from sympy.utilities.pytest import XFAIL, skip, raises

n= Symbol('n')
a = Symbol('a')
b = Symbol('b')
C0 = Symbol('C0')
C1 = Symbol('C1')
C2 = Symbol('C2')

f = Function('f')
g = Function('g')
K = Function('K')

vfnh = ['Volterra', 'First Kind', 'Non-homogenous']
vfh = ['Volterra', 'First Kind', 'Homogenous']
vsnh = ['Volterra', 'Second Kind', 'Non-homogenous']
vsh = ['Volterra', 'Second Kind', 'Homogenous']
ffnh = ['Fredholm', 'First Kind', 'Non-homogenous']
ffh = ['Fredholm', 'First Kind', 'Homogenous']
fsnh = ['Fredholm', 'Second Kind', 'Non-homogenous']
fsh = ['Fredholm', 'Second Kind', 'Homogenous']

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
    
def test_solveapproximate():
    eq1 = Eq((5.0/6.0)*x+Integral(0.5*x*y*f(y),(y,0,1)),f(x)) # Example in maple's intsolve
    assert solve_approximate(eq1,f(x),21) == x
    eq2 = Eq(1 + Integral(f(y),(y,0,x)),f(x))
    assert solve_approximate(eq2,f(x),5) == 1 + x + x**2/2\
     + x**3/6 + x**4/24 + x**5/120 # Series expansion for exp(x)
    eq3 = Eq(1 + x + Integral(2*x*f(y),(y,0,1)),f(x))
    assert solve_approximate(eq3,f(x),2,x) == 1 + 5*x
    eq4 =  Eq(1 + Integral(x*x*y*f(y),(y,0,1)),f(x))
    assert solve_approximate(eq4,f(x),1,x) == 1 + x**2/3
    eq5 =  Eq(1 + x + Integral(-f(y),(y,0,x)),f(x)) # Converges to 1
    assert solve_approximate(eq5,f(x),1,x) == 1 + x**7/5040
    eq6 = Eq(1 + Integral(x*f(y),(y,0.0,1.0)),f(x)) 
    assert solve_approximate(eq6,f(x),49) == 1 + 2.0*x
    
def test_solveadomian():
    eq1 = Eq(1 + Integral(f(y),(y,0,x)),f(x))
    assert solve_approximate(eq1,f(x),6) == 1 + x + x**2/2\
     + x**3/6 + x**4/24 + x**5/120 # Series expansion for exp(x)
    eq2 =  Eq(exp(x) + x*exp(x) + Integral(-exp(x-y)*f(y),\
                                           (y,0,x)),f(x))
    eq3 = Eq(1 + x + Integral(-f(y),(y,0,x)),f(x))
    assert solve_adomian(eq1, f(x), 9) == 1 + x**9/362880 # Converges to 1

def test_solveseries():
    eq1 =  Eq(1 + Integral(f(y),(y,0,1)),f(x))
    assert solve_series(eq1,f(x),3) == C0 - 3*x**2
    eq2 =  Eq(1 + Integral((x**2-y**2)*f(y),(y,0,1)),f(x))
    assert solve_series(eq1,f(x),10) == 30/49 + 45*x**2/49
    eq3 =  Eq(1 + Integral(f(y),(y,0,x)),f(x))
    assert solve_series(eq3,f(x),10) == 1 + x + x**2/2 + \
                                        x**3/6 + x**4/24 + x**5/120
def test_checkidesol():
    eq1 =  Eq(1 + Integral(f(y),(y,0,x)),f(x))
    assert checkidesol(eq1,f(x),exp(x))     
    eq2 = Eq(1 + Integral(x*f(y),(y,0.0,1.0)),f(x))
    assert checkidesol(eq1, f(x),solve_approximate(eq1,f(x),49))
