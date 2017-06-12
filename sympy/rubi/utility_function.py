'''
Utility functions for Constraints in Rubi
'''

from sympy.functions import (log, sin, cos, sqrt)
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify
from sympy.functions.elementary.integers import floor, frac
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acsch
from sympy.functions.elementary.trigonometric import atan, acsc, asin, asin, acos
from sympy.polys.polytools import degree, Poly
from sympy.simplify.simplify import fraction, simplify, count_ops
from sympy.integrals.integrals import integrate
from fractions import Fraction
from mpmath import hyp2f1, ellippi, ellipe, ellipf, appellf1

def ZeroQ(expr):
    return expr == 0

def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    if isinstance(nodes, list):
        return not any(expr.has(var) for expr in nodes)
    else:
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
    else:
        return False

def Sqrt(a):
    return sqrt(a)

def ArcCosh(a):
    return acosh(a)

def Coefficient(expr, var, n):
    a = Poly(expr, var)
    if (degree(a) - n) < 0:
        return 0
    else:
        return a.all_coeffs()[degree(a) - n]

def RemoveContent(expr, x):
    for i in expr.args:
        if not i.has(x):
            expr = expr - i
    return expr

def Denominator(var):
    try:
        return fraction(sympify(Fraction (var)), exact=True)[1]
    except TypeError:
        return fraction(var, exact=True)[1]

def Hypergeometric2F1(a, b, c, z):
    return hyp2f1(a, b, c, z)

def ArcTan(a):
    return atan(a)

def Not(var):
    return not(var)

def Simplify(expr):
    return simplify(expr)

def FractionalPart(a):
    return FracPart(a)

def IntegerPart(a):
    return IntPart(a)

def SumSimplerQ(u, v):
    # If u+v is simpler than u, SumSimplerQ(u,v) returns True, else it returns False
    return count_ops(u+v) < count_ops(u)

def SimplerQ(u, v):
    # If u is simpler than v, SimplerQ(u,v) returns True, else it returns False. SimplerQ(u,u) returns False.
    return count_ops(u) < count_ops(v)

# utility functions used in RUBI tests

def AppellF1(a, b1, b2, c, x, y):
    return appellf1(a, b1, b2, c, x, y)

def Integrate(f, x):
    return integrate(f, x)

def EllipticPi(*args):
    return ellippi(*args)

def EllipticE(*args):
    return ellipe(*args)

def EllipticF(Phi, m):
    return ellipf(Phi, m)

def ArcTan(a):
    return atan(a)

def ArcTanh(a):
    return atanh(a)

def ArcSin(a):
    return asin(a)

def ArcSinh(a):
    return asinh(a)

def ArcCos(a):
    return acos(a)

def ArcCsc(a):
    return acsc(a)

def ArcCsch(a):
    return acsch(a)
