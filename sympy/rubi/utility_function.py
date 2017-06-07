'''
Utility functions for Constraints in Rubi
'''

from sympy.functions import (log, sin, cos, sqrt)
from sympy.core.symbol import Symbol
from sympy.functions.elementary.integers import floor, frac
from sympy.functions.elementary.hyperbolic import acosh
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
