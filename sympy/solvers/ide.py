"""
This module contains iesolve() and different helper functions that it
uses.

iesolve() solves linear integral equations. It replicates the functionality of intsolve in Maple.

**Functions in this module**

    These are the user functions in this module:

    - iesolve() - Solves IDEs. It takes the integral equation and the unknown function as its parameters
                  as well as the method to be used to solve the integral equation. This method only solves
                  linear integral equations.
    - iesolve_nonlinear() - Solves certain kinds of non linear integral equations.              
    - classify_ide() - Classifies IDEs into possible hints for iesolve(). 
    - checkidesol() - Checks if the given function is the solution of an IDE.
    
    References
    [1] Integral Equations and their Applications, M. Rahman
    [2] Handbook of Integral Equations, Andrei D. Polyanin and Alexander V. Manzhirov
"""

from sympy.core import Add, Basic, C, S, Mul, Pow, oo
from sympy.core.function import Derivative, diff, expand_mul, FunctionClass
from sympy.core.relational import Equality, Eq
from sympy import Function, simplify
from sympy.functions import cos, exp, im, log, re, sin, sign

from sympy.solvers import solve, dsolve
from sympy import (S, symbols, Integral, Derivative, integrate, exp, oo, Symbol,
        Function, Rational, log, sin, cos, pi, E, I, Poly, LambertW, diff,
        sympify, sqrt, atan, asin, acos, atan, DiracDelta, Heaviside,
        Lambda)

from sympy.abc import x, y, z

def classify_ide_type(term, classifier, symbol):
    """
    Classifies an integral equation as Fredholm or Volterra by
    checking if integral has definite or indefinite limits
    """
    eqlimits = term.limits[0][1]
    if symbol in eqlimits:
        classifier.insert(0,'Volterra')
    else:
        classifier.insert(0,'Fredholm')
    return classifier

def classify_ide_kind(term, classifier, func):
    #Resolves the kind of integral equation
    classifier = classify_ide_type(term, classifier, func.args[0])
    if term.function.func == func.func:
        classifier.insert(1,'First Kind')
    else:
        for mult in term.function.args:
            # We compare function classes and not the functions themselves
            if mult.func == func.func:
                classifier.insert(1,'First Kind')
    
    if not ('First Kind') in classifier:
        raise ValueError("The unknown function does not appear in the integral")              
    return classifier

def process_terms(term, func, classifier):
    """
    A recursive solution is required to take care of nested Muls
    in an expression
    """
    if type(term) == Integral:
        return classify_ide_kind(term, classifier, func)
    
    if term == func or term.func == func:
        classifier.insert(1,'Second Kind')
        return classifier
    else:
        if not ('Non-homogenous') in classifier:
            classifier.insert(2,'Non-homogenous')

    if type(term) == Mul:
        for mult in term.args: 
            classifier = process_terms(mult, func, classifier)
                            
    return classifier         
        
def classify_ide(eq, func, dicr=False):
    """
    Returns a list indicating the classification of the given integral equation
    
    **Examples**
    
        >>> from sympy import Function, classify_ode, Eq
        >>> from sympy.abc import x
    
        Consider the following integral equation
        >>> eq = Eq(Integral(f(y)*exp(x-y),(y,0,x)),x)
        >>> classify_ide (eq, f(x))
        ['Volterra', 'First Kind', 'Non-homogenous']
        
        Define an arbitrary function K
        >>> K = Function('K')
        >>> eq = Eq(Integral(f(y)*K(x-y),(y,0,x)),f(x))
        >>> classify_ide (eq, f(x))
        ['Volterra', 'Second Kind', 'Homogenous']
        
        >>> eq =  Eq(f(x) + n*Integral(f(y)*((x-y)**(0.5)),(y,0,x)),g(x))
        >>> classify_ide (eq, f(x))
        ['Non-homogenous', 'Second Kind', 'Volterra']
    """
    if len(func.args) != 1:
        raise ValueError("iesolve() and classify_ide() only work with functions " + \
            "of one variable")
    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return classify_ide(eq.lhs-eq.rhs, func)
        eq = eq.lhs
    
    classifier = [""]
    if type(eq) == Integral:
        # The corner case of an integral equation of the first kind where known function is null
        classifier.insert(1,'First Kind')
        classifier.insert(0,'Fredholm')
        # Check if Volterra or Fredholm
        classifier = classify_ide_type(eq, classifier, func.args[0])
    else:                
        for term in eq.args:
            process_terms(term, func, classifier)    
        if ('Second Kind') in classifier:
            classifier.remove('First Kind')
            if len(eq.args) < 3:
                classifier.remove('Non-homogenous')
                classifier.insert(2,'Homogenous')    
    
    classifier.remove('')
    return classifier

def checkidesol(ide, func, fn):
    """
    This method verifies if the given function satisfies the given Integral equation
    """
    eq = ide
    if len(func.args) != 1:
        raise ValueError("checkidesol() will only work with functions " + \
            "of one variable")
    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return checkidesol(eq.lhs-eq.rhs, func, fn)
        eq = eq.lhs
        
    neweq = Eq(0,0)
    
    #We need to do termwise replacement here. 
    #This code will be updated when polysys12 is pushed into the main repo
    for term in eq.args:
        if type(term) == Integral:
            integralsymbol = term.variables[0]
            subsfunc = fn.subs(func.args[0], integralsymbol)
            neweq = Eq(neweq.lhs + integrate(term.function.subs(f(y), subsfunc), term.limits))
        else:
            neweq = Eq(neweq.lhs + term.subs(func, fn),0)
    return neweq == Eq(0,0)        

def iesolve(eq, func, method = ""):
    """
    This is the main routine exposed to the user for solving linear integral
    equations. This routine first attempts to classify
    the given integral equation. It then selects the method most appropriate for
    solving the integral equation.
    """
    ide_classification = classify_ide(eq, func)
    raise NotImplementedError()

def iesolve_nonlinear(eq, func, method = ""):
    raise NotImplementedError()

def subs_func_in_integral(func, term, neweq, startsoln):
    # Helper function to substitute the starting solution in the integral term
    funcsymbol = func.args[0] #x
    integralsymbol = term.variables[0] #y
    subsfunc = func.subs(funcsymbol, integralsymbol) #f(y)
    try:
        startsolnsubs = startsoln.subs(funcsymbol, integralsymbol)
    except AttributeError:
        startsolnsubs = startsoln
    neweq = Eq(neweq.lhs + integrate(term.function.subs(subsfunc, startsolnsubs),\
                                      term.limits),neweq.rhs)
    return neweq

def solve_adomian(eq, func, n):
    """
    Solves integral equations using the Adomian Decomposition technique
    The user needs to specify the number of terms upto which the series
    needs to be calculated. The non integral terms are isolated and these are then
    assigned to U_0(x). The rest of the terms are calculated using the recursive
    relation U_n(x) = Integral(K(x,t)*U_(n-1)(t), (t,0,x)), where K(x,t) is a known
    function within the integral term. K is also known as the kernel of the equation
    """
    components = []
    for i in range(n):
        component_name = "U_" + str(i)
        components.insert(i,(Function(component_name))(func.args[0]))

    series = Eq(0,0)
    for i in range(n):
        series = Eq(series.lhs + components[i])
    series = series.lhs
    
    neweq = Eq(0,0)
    nonintegral = Eq(0,0)
    if type(eq.lhs) == Integral:
        termsymbol = eq.lhs.variables[0]
        neweq = subs_func_in_integral(func, eq.lhs, neweq, series)
        subsfunc = func.subs(func.args[0], termsymbol)
        kernel = eq.lhs.function.subs(subsfunc, 1)
    else:    
        for term in eq.lhs.args:
            if type(term) == Integral:
                termsymbol = term.variables[0]
                seriesfunc = series.subs(func.args[0], termsymbol)
                subsfunc = func.subs(func.args[0], termsymbol)
                kernel = term.function.subs(subsfunc, 1)
                neweq = Eq(neweq.lhs + Integral(term.function.subs(subsfunc, seriesfunc),\
                                      term.limits),neweq.rhs)
            else:
                neweq = Eq(neweq.lhs + term.subs(func, series), neweq.rhs)
                nonintegral = Eq(nonintegral.lhs + term)
    # Isolate U_0(x)
    components[0] = nonintegral.lhs
    
    # Calculate the rest of the functions
    for i in range(1,n):
        if type(components[i-1]) != int:
            components[i] = integrate(kernel*components[i-1].subs(func.args[0],termsymbol),(termsymbol,0,func.args[0]))
        else:    
            components[i] = integrate(kernel*components[i-1],(t,0,func.args[0]))
    ans = Eq(0,0)
    for i in range(n):
        ans = Eq(ans.lhs + components[i])    
    return ans.lhs

def solve_series(eq, func, n):
    """
    This method assumes that the solution we are going to get
    is analytic and hence can be expressed as a power series. This is then
    substituted in the equation for the unknown function and evaluated
    term by term. Coefficients of powers of x are collected on one side
    and then a series of linear equations are obtained which are solved
    to give the values of the coefficients in the power series. 
    Use for Volterra integral equations of the second kind where K(x,x) = 0.
    """
    components = []
    symbols = []
    for i in range(0,n,1):
        component_name = "U_" + str(i)
        term_coeff = Symbol("C" + str(i))
        symbols.insert(i, term_coeff)
        components.insert(i, (Function(component_name))(func.args[0]))
        ith_term = term_coeff*func.args[0]**i
        components[i] = Eq(components[i], ith_term)
        
    series = Eq(0,0)
    for i in range(0,n,1):
        series = Eq(series.lhs + components[i].rhs)
    series = series.lhs
    
    neweq = Eq(0,eq.rhs.subs(func, series))
    nonintegral = Eq(0,0)
    for term in eq.lhs.args:
        if type(term) == Integral:
            termsymbol = term.variables[0]
            seriesfunc = series.subs(func.args[0], termsymbol)
            subsfunc = func.subs(func.args[0], termsymbol)
            neweq = Eq(neweq.lhs + (integrate(term.function.subs(subsfunc, seriesfunc),\
                                  term.limits)).expand(),neweq.rhs) 
        else:
            nonintegral = Eq(nonintegral.lhs + term.expand(),0)
    
    neweq = Eq(collect(neweq.lhs + nonintegral.lhs - neweq.rhs,\
                       [func.args[0]**i for i in range(1,n,1)]),0)

    series_system = []
    series_system.insert(0, Eq(0))
    for i, term in zip(range(len(neweq.lhs.args)), neweq.lhs.args):        
        if type(term) == Mul and func.args[0] in term.args[0]:
            series_system.append(Eq(term.args[1]))
        else:
            if func.args[0] in term:
                if term.args[1] == symbols[len(symbols)-1]:
                    continue
                series_system.append(Eq(term.subs(func.args[0],1)))
                continue
            series_system[0] = Eq(series_system[0].lhs + term)
    
    if len(series_system) > len(symbols):
        del series_system[len(symbols):len(series_system)]    
    
    series_soln = solve(series_system, symbols)
    if series_soln == None:
        raise ValueError ("No power series solution exists")
    for k, v in series_soln.iteritems():
        series = series.subs(k, v)
    
    return series

def solve_approximate(eq, func, level, startsoln = 1, picard = False):
    """
    Gives an approximate solution to an integral equation. We plugin the starting
    solution into the integral term and solve. This solution is then used for the 
    next iteration. Level controls the number of iterations. 
    It is expected that the integral equation is given in an explicit form, with
    the RHS having func. Valid only for second kind integral equations.
    """
    if picard:
        nonintegral = Eq(0,0)
        for term in eq.lhs.args:
            if type(term) != Integral:
                nonintegral = Eq(nonintegral.lhs + term)
        startsoln = nonintegral
                            
    if level < 1:
        return startsoln
    else:
        neweq = Eq(0,eq.rhs)
        if type(eq.lhs) == Integral:
            neweq = subs_func_in_integral(func, eq.lhs, neweq, startsoln)
        else:    
            for term in eq.lhs.args:
                if type(term) == Integral:
                    neweq = subs_func_in_integral(func, term, neweq, startsoln)
                else:
                    neweq = Eq(neweq.lhs + term.subs(func, startsoln), neweq.rhs)        
    return solve_approximate(eq, func, level - 1, neweq.lhs, False)        

def solve_asode(eq, func, maxdepth, depth = 1, initialvalues=[]):
    """
    Differentiate the equation to get a differential equation and use dsolve
    It is impractical and at other times impossible to solve certain IDEs by reducing
    them to ODEs. To counter this situation we enforce a depth constraint which restricts
    the number of times an integral equation can be differentiated to get a pure
    differential equation. This method is only valid for volterra integral equations
    
    *Warning*
    The variable maxdepth should be set to a reasonably low value like 5-6 else
    the routine can take an inordinately large amount of time to finish.
    
    """
    #FIXME: This gives incorrect answers because we dont take into account initial value
    #Also the initial values are at the lower limit of the integral term and not always at 0
    if depth > maxdepth:
        raise ValueError ("Solving this integral equation as an ODE is not feasible")
    
    neweq = Eq(0,eq.rhs)
    needs_diff = False
    for term in eq.lhs.args:
        if type(term) == Integral:
            term_lowerlt = term.limits[0][1][0]
        diffterm = diff(term, func.args[0],1)
        if type(diffterm) == Integral:
            term_lowerlt = diffterm.limits[0][1][0]
            needs_diff = True
        if type(diffterm) == Mul or type(diffterm) == Add:
            for innerterm in diffterm.args:
                if type(innerterm) == Integral:
                    term_lowerlt = innerterm.limits[0][1][0]
                    needs_diff = True            
        neweq = Eq(neweq.lhs + diffterm,neweq.rhs)
    neweq = Eq(neweq.lhs, diff(neweq.rhs, func.args[0],1))
    if len(initialvalues) == 0:
        initialvalues.append(Eq(eq.lhs.subs(func.args[0], term_lowerlt),\
                                eq.rhs.subs(func.args[0], term_lowerlt)))
    initialvalues.append(Eq(neweq.lhs.subs(func.args[0], term_lowerlt),\
                            neweq.rhs.subs(func.args[0], term_lowerlt)))
  
    if needs_diff:
        neweq = solve_asode(neweq, func, depth + 1, initialvalues)
        return neweq
    #NOTE: We need to make use of the initialvalues somehow
    return dsolve(neweq, func)

def solve_neumann(eq,func,n):
    """
    Solves integral equations using Neumann Series method also 
    known as the method of successive substitutions. 
    """
    #Create terms of Neumann series
    terms = []
    for i in range(0, n):
        terms.insert(i, Eq(Function("U" + str(i))(func.args[0]),0))
    non_integral_terms = Eq(0,0)
    
    #Extract kernel and non-integral terms
    for term in eq.lhs.args:
        if type(term) != Integral:
            non_integral_terms = Eq(non_integral_terms.lhs + term, 0)
        if type(term) == Integral:
            limits = term.limits
            termsymbol = term.variables[0]
            subsfunc = func.subs(func.args[0], termsymbol)
            kernel = term.function.subs(subsfunc, 1)
    
    #For each term find out the required integral and print the last term
    terms[0] = Eq(terms[0].lhs,\
                   integrate(kernel * non_integral_terms.lhs, limits))
    
    for i in range(1, n):
        terms[i] = Eq(terms[i].lhs, \
                      integrate(kernel * terms[i-1].rhs.subs(func.args[0], termsymbol), limits))
                        
    ans = Eq(0,0)
    for i in range(n):
        ans = Eq(ans.lhs + terms[i].rhs)    
    return ans.lhs 
    
def solve_eigenfunction():
    """
    Solves integral equations using the Eigenfunction technique
    """
    raise NotImplementedError()

def solve_laplace():
    """
    Solves integral equations using Laplace Transform technique
    """
    raise NotImplementedError()
