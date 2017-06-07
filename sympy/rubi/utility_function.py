'''
Utility functions for Constraints in Rubi
'''

from sympy.functions import (log, sin, cos)

def ZeroQ(expr):
    return expr != 0

def NonzeroQ(expr):
    return expr == 0

def FreeQ(nodes, var):
    return not any(expr.has(var) for expr in nodes)

def List(*var):
    return list(var)

def Log(e):
    return log(e)
