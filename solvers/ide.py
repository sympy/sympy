"""
This module contains iesolve() and different helper functions that it
uses.

iesolve() solves linear integral equations. It replicates the functionality of intsolve in Maple.

**Functions in this module**

    These are the user functions in this module:

    - iesolve() - Solves IDEs. It takes the integral equation and the unknown function as its parameters
                  as well as the method to be used to solve the integral equation.
    - classify_ide() - Classifies IDEs into possible hints for iesolve(). 
    - checkidesol() - Checks if the given function is the solution of an IDE.
    
"""

from sympy.core import Add, Basic, C, S, Mul, Pow, oo
from sympy.core.function import Derivative, diff, expand_mul, FunctionClass
from sympy.core.multidimensional import vectorize
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Symbol, Wild
from sympy.core.sympify import sympify

from sympy import Function
from sympy.functions import cos, exp, im, log, re, sin, sign
from sympy.simplify import collect, logcombine, powsimp, separatevars, \
    simplify, trigsimp
from sympy.solvers import solve
from sympy import (S, symbols, integrate, Integral, Derivative, exp, oo, Symbol,
        Function, Rational, log, sin, cos, pi, E, I, Poly, LambertW, diff,
        Matrix, sympify, sqrt, atan, asin, acos, atan, DiracDelta, Heaviside,
        raises, Lambda, sstr)

from sympy.utilities import numbered_symbols, all, any, make_list
from sympy.abc import x, y, z
from collections import Hashable

f = Function('f')       
sol1 = Eq(Integral(f(y)*exp(x-y),(y,0,x)),x)

def classify_ide_type(term, classifier):
    eqlimits = term.limits[0][1]
    if eqlimits.__contains__(term.limits[0][0]):
        classifier.append('Volterra')
    else:
        classifier.append('Fredholm')
    return classifier
        
def classify_ide(eq, func, dicr=False):
    """
    Returns a list indicating the classification of the given integral equation
    """
    if len(func.args) != 1:
        raise ValueError("iesolve() and classify_ide() only work with functions " + \
            "of one variable")
    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return classify_ide(eq.lhs-eq.rhs, func)
        eq = eq.lhs
    
    classifier = []
    if type(eq) == Integral:
        # The corner case of an integral equation of the first kind where known function is null
        classifier.append('First Kind')
        classifier.append('Fredholm')
        # Check if Volterra or Fredholm
        classifier = classify_ide_type(eq, classifier)
    else:                
        for term in eq.args:
            if type(term) == Integral:
                # Check if Volterra or Fredholm
                classifier = classify_ide_type(term, classifier)

                if term.function == func:
                    classifier.append('First Kind')
                else:
                    for mult in term.function.args:
                        # We compare function classes and not the functions themselves
                        if mult.func == func.func:
                            classifier.append('First Kind')
                
                if not ('First Kind') in classifier:
                    raise ValueError("The unknown function does not appear in the integral")              
                continue
            
            if term == func or term.func == func:
                classifier.append('Second Kind')
            else:
                if not ('Non-homogenous') in classifier:
                    classifier.append('Non-homogenous')            
                     
        if ('Second Kind') in classifier:
            classifier.remove('First Kind')
            if len(eq.args) < 3:
                classifier.append('Homogenous')    
    
    return classifier

def check_idesol(eq, func, __f__):
    """
    This method verifies if the given function satisfies the given Integral equation
    """
    if len(func.args) != 1:
        raise ValueError("check_idesol() will only work with functions " + \
            "of one variable")
    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return check_idesol(eq.lhs-eq.rhs, func, __f__)
        
    if type(__f__.func) == FunctionClass:
        __f__ = __f__.func
    neweq = Eq(0,0)
    """
    We need to do termwise replacement here. 
    This code will be updated when polysys12 is pushed into the main repo
    """
    for term in eq.args:
        if type(term) == Integral:
            integralsymbol = term.variables[0]
            subsfunc = func.subs(func.args[0], term.variables[0]) #f(y)
            neweq = Eq(neweq.lhs + integrate(term.function.subs(subsfunc, __f__.subs(func.args[0], integralsymbol)), term.limits),0)
        else:
            neweq = Eq(neweq.lhs + term.subs(func, __f__),0)
    return neweq == Eq(0,0)        
