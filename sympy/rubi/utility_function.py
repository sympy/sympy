'''
Utility functions for Constraints in Rubi
'''

from sympy.functions import (log, sin, cos, sqrt)
from sympy.core.symbol import Symbol
from sympy.functions.elementary.integers import floor, frac
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acsch
from sympy.functions.elementary.trigonometric import atan, acsc, asin, asin
from sympy.polys.polytools import degree

def ZeroQ(expr):
    return expr == 0

def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    if isinstance(nodes, list):
        return not any(expr.has(var) for expr in nodes)
    elif isinstance(nodes, Symbol):
        return nodes != var

def List(*var):
    return list(var)

def Log(e):
    return log(e)

def PositiveIntegerQ(var):
    return var.is_Integer and var > 0

def NegativeIntegerQ(var):
    return var.is_Integer and var < 0

def PositiveQ(var):
    return var > 0

def IntegerQ(var):
    return var.is_Integer

def PosQ(var):
    return var > 0

def NegQ(var):
    return var < 0

def FracPart(var):
    return frac(var)

def IntPart(var):
    return floor(var)

def RationalQ(var):
    return var.is_Rational

def Subst(a, x, y):
    return a.subs(x, y)

def LinearQ(expr, x):
    if degree(expr, gen=x) == 1:
        return True

def Sqrt(a):
    return sqrt(a)

def ArcCosh(a):
    return acosh(a)

def TogetherSimplify(expr):
    return

def Coefficient(u, var, n):
    return

def RemoveContent():
    return

def Sqrt(a):
    return math.sqrt(a)

def ExpandIntegrand(expr, x):
    return

def ArcCosh(a):
    return acosh(a)

def With():
    return

def Denominator(var):
    return fraction(var, exact=True)[1]

def Hypergeometric2F1(a, b, c, z):
    return hyp2f1(a, b, c, z)

def TogetherSimplify():
    return

def IntLinearcQ():
    return

def ArcTan(a):
    return atan(a)

def Not(var):
    return not(var)

def Simplify():
    return

def FractionalPart(a):
    return FracPart(a)

def IntegerPart(a):
    return IntPart(a)

def Simp():
    return

def Rt():
    return

def SumSimplerQ():
    return

# utility functions used in tests

def AppellF1(a, b1, b2, c, x, y):
    return appellF1(a, b1, b2, c, x, y)

def Integrate(f, x):
    return integrate(f, x)

def hypergeom():
    return

def EllipticPi(*args):
    return ellippi(*args)

def EllipticE(*args):
    return ellipe(*args)

def EllipticF(Phi, m):
    return ellipf(Phi, m)

def arctan(a):
    return atan(a)

def arctanh(a):
    return atanh(a)

def arcsin(a):
    return asin(a)

def arcsinh(a):
    return asinh(a)

def arccos(a):
    return acos(a)

def arccosh(a):
    return acosh(a)

def arccsc(a):
    return acsc(a)

def arccsch(a):
    return acsch(a)
