from sympy import S, symbols
from sympy.stats import Poisson, Geometeric
from sympy.stats.joint_rv import Joint

x, y = symbols('x y')

def test_density():
    from sympy.functions.elementary.exponential import exp
    from sympy.functions.combinatorial.factorials import factorial
    X, Y = Geometric('X', S(1)/2), Poisson('Y', 4)
    Z = Joint('Z', (X, Y))
    assert density(Z) == 2**(-x + 1)*4**y*exp(-4)/(2*factorial(y))