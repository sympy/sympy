""" This module contain solvers for all kinds of equations:

    - algebraic, use solve()

    - recurrence, use rsolve()

    - differential, use dsolve()

    - transcendental, use tsolve()

    - nonlinear (numerically), use nsolve()
      (you will need a good starting point)

"""

from sympy.core.sympify import sympify
from sympy.core.basic import Basic, S, C, Mul
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.core.symbol import Symbol, Wild
from sympy.core.relational import Equality
from sympy.core.function import Derivative, diff, Function, expand_mul
from sympy.core.numbers import ilcm
from sympy.core.multidimensional import vectorize

from sympy.functions import sqrt, log, exp, LambertW, sin, cos, re, im
from sympy.simplify import simplify, collect, logcombine, separatevars, trigsimp
from sympy.matrices import Matrix, zeros, wronskian
from sympy.polys import roots, RootsOf, discriminant, RootOf

from sympy.utilities import any, all, numbered_symbols
from sympy.utilities.lambdify import lambdify
from sympy.mpmath import findroot

from sympy.solvers.polysys import solve_poly_system

from warnings import warn

# Codes for guess solve strategy
GS_POLY = 0
GS_RATIONAL = 1
GS_POLY_CV_1 = 2 # can be converted to a polynomial equation via the change of variable y -> x**a, a real
GS_POLY_CV_2 = 3 # can be converted to a polynomial equation multiplying on both sides by x**m
                 # for example, x + 1/x == 0. Multiplying by x yields x**2 + x == 0
GS_RATIONAL_CV_1 = 4 # can be converted to a rational equation via the change of variable y -> x**n
GS_TRANSCENDENTAL = 5

def guess_solve_strategy(expr, symbol):
    """
    Tries to guess what approach should be used to solve a specific equation

    Returns
    =======
       - -1: could not guess
       - integer > 0: code representing certain type of equation. See GS_* fields
         on this module for a complete list

    Examples
    ========
    >>> from sympy import Symbol, Rational
    >>> x = Symbol('x')
    >>> guess_solve_strategy(x**2 + 1, x)
    0
    >>> guess_solve_strategy(x**Rational(1,2) + 1, x)
    2
    """
    eq_type = -1
    if expr.is_Add:
        return max([guess_solve_strategy(i, symbol) for i in expr.args])

    elif expr.is_Mul:
        # check for rational functions
        num, denom = expr.as_numer_denom()
        if denom != 1 and denom.has(symbol):
            #we have a quotient
            m = max(guess_solve_strategy(num, symbol), guess_solve_strategy(denom, symbol))
            if m == GS_POLY:
                return GS_RATIONAL
            elif m == GS_POLY_CV_1:
                return GS_RATIONAL_CV_1
            else:
                raise NotImplementedError
        else:
            return max([guess_solve_strategy(i, symbol) for i in expr.args])

    elif expr.is_Symbol:
        return GS_POLY

    elif expr.is_Pow:
        if expr.exp.has(symbol):
            return GS_TRANSCENDENTAL
        elif not expr.exp.has(symbol) and expr.base.has(symbol):
            if expr.exp.is_Integer and expr.exp > 0:
                eq_type = max(eq_type, GS_POLY)
            elif expr.exp.is_Integer and expr.exp < 0:
                eq_type = max(eq_type, GS_POLY_CV_2)
            elif expr.exp.is_Rational:
                eq_type = max(eq_type, GS_POLY_CV_1)
            else:
                return GS_TRANSCENDENTAL

    elif expr.is_Function and expr.has(symbol):
        return GS_TRANSCENDENTAL

    elif not expr.has(symbol):
        return GS_POLY

    return eq_type

def solve(f, *symbols, **flags):
    """Solves equations and systems of equations.

       Currently supported are univariate polynomial and transcendental
       equations and systems of linear and polynomial equations.  Input
       is formed as a single expression or an equation,  or an iterable
       container in case of an equation system.  The type of output may
       vary and depends heavily on the input. For more details refer to
       more problem specific functions.

       By default all solutions are simplified to make the output more
       readable. If this is not the expected behavior,  eg. because of
       speed issues, set simplified=False in function arguments.

       To solve equations and systems of equations of other kind, eg.
       recurrence relations of differential equations use rsolve() or
       dsolve() functions respectively.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       Solve a polynomial equation:

       >>> solve(x**4-1, x)
       [1, -1, -I, I]

       Solve a linear system:

       >>> solve((x+5*y-2, -3*x+6*y-15), x, y)
       {x: -3, y: 1}

    """
    if not symbols:
        raise ValueError('no symbols were given')

    if len(symbols) == 1:
        if isinstance(symbols[0], (list, tuple, set)):
            symbols = symbols[0]

    symbols = map(sympify, symbols)
    result = list()

    # Begin code handling for Function and Derivative instances
    # Basic idea:  store all the passed symbols in symbols_passed, check to see
    # if any of them are Function or Derivative types, if so, use a dummy
    # symbol in their place, and set symbol_swapped = True so that other parts
    # of the code can be aware of the swap.  Once all swapping is done, the
    # continue on with regular solving as usual, and swap back at the end of
    # the routine, so that whatever was passed in symbols is what is returned.
    symbols_new = []
    symbol_swapped = False

    if isinstance(symbols, (list, tuple)):
        symbols_passed = symbols[:]
    elif isinstance(symbols, set):
        symbols_passed = list(symbols)

    i = 0
    for s in symbols:
        if s.is_Symbol:
            s_new = s
        elif s.is_Function:
            symbol_swapped = True
            s_new = Symbol('F%d' % i, dummy=True)
        elif s.is_Derivative:
            symbol_swapped = True
            s_new = Symbol('D%d' % i, dummy=True)
        else:
            raise TypeError('not a Symbol or a Function')
        symbols_new.append(s_new)
        i += 1

        if symbol_swapped:
            swap_back_dict = dict(zip(symbols_new, symbols))
    # End code for handling of Function and Derivative instances

    if not isinstance(f, (tuple, list, set)):
        f = sympify(f)

        # Create a swap dictionary for storing the passed symbols to be solved
        # for, so that they may be swapped back.
        if symbol_swapped:
            swap_dict = zip(symbols, symbols_new)
            f = f.subs(swap_dict)
            symbols = symbols_new

        if isinstance(f, Equality):
            f = f.lhs - f.rhs

        if len(symbols) != 1:
            raise NotImplementedError('multivariate equation')

        symbol = symbols[0]

        strategy = guess_solve_strategy(f, symbol)

        if strategy == GS_POLY:
            poly = f.as_poly( symbol )
            assert poly is not None
            result = roots(poly, cubics=True, quartics=True).keys()

        elif strategy == GS_RATIONAL:
            P, Q = f.as_numer_denom()
            #TODO: check for Q != 0
            result = solve(P, symbol, **flags)

        elif strategy == GS_POLY_CV_1:
            args = list(f.args)
            if isinstance(f, Add):
                # we must search for a suitable change of variable
                # collect exponents
                exponents_denom = list()
                for arg in args:
                    if isinstance(arg, Pow):
                        exponents_denom.append(arg.exp.q)
                    elif isinstance(arg, Mul):
                        for mul_arg in arg.args:
                            if isinstance(mul_arg, Pow):
                                exponents_denom.append(mul_arg.exp.q)
                assert len(exponents_denom) > 0
                if len(exponents_denom) == 1:
                    m = exponents_denom[0]
                else:
                    # get the GCD of the denominators
                    m = ilcm(*exponents_denom)
                # x -> y**m.
                # we assume positive for simplification purposes
                t = Symbol('t', positive=True, dummy=True)
                f_ = f.subs(symbol, t**m)
                if guess_solve_strategy(f_, t) != GS_POLY:
                    raise TypeError("Could not convert to a polynomial equation: %s" % f_)
                cv_sols = solve(f_, t)
                for sol in cv_sols:
                    result.append(sol**m)

            elif isinstance(f, Mul):
                for mul_arg in args:
                    result.extend(solve(mul_arg, symbol))

        elif strategy == GS_POLY_CV_2:
            m = 0
            args = list(f.args)
            if isinstance(f, Add):
                for arg in args:
                    if isinstance(arg, Pow):
                        m = min(m, arg.exp)
                    elif isinstance(arg, Mul):
                        for mul_arg in arg.args:
                            if isinstance(mul_arg, Pow):
                                m = min(m, mul_arg.exp)
            elif isinstance(f, Mul):
                for mul_arg in args:
                    if isinstance(mul_arg, Pow):
                        m = min(m, mul_arg.exp)
            f1 = simplify(f*symbol**(-m))
            result = solve(f1, symbol)
            # TODO: we might have introduced unwanted solutions
            # when multiplied by x**-m

        elif strategy == GS_TRANSCENDENTAL:
            #a, b = f.as_numer_denom()
            # Let's throw away the denominator for now. When we have robust
            # assumptions, it should be checked, that for the solution,
            # b!=0.
            result = tsolve(f, *symbols)
        elif strategy == -1:
            raise ValueError('Could not parse expression %s' % f)
        else:
            raise NotImplementedError("No algorithms are implemented to solve equation %s" % f)

        if symbol_swapped:
            result = [ri.subs(swap_back_dict) for ri in result]

        if flags.get('simplified', True) and strategy != GS_RATIONAL:
            return map(simplify, result)
        else:
            return result
    else:
        if not f:
            return {}
        else:
            # Create a swap dictionary for storing the passed symbols to be
            # solved for, so that they may be swapped back.
            if symbol_swapped:
                swap_dict = zip(symbols, symbols_new)
                f = [fi.subs(swap_dict) for fi in f]
                symbols = symbols_new

            polys = []

            for g in f:
                g = sympify(g)

                if isinstance(g, Equality):
                    g = g.lhs - g.rhs

                poly = g.as_poly(*symbols)

                if poly is not None:
                    polys.append(poly)
                else:
                    raise NotImplementedError()

            if all(p.is_linear for p in polys):
                n, m = len(f), len(symbols)
                matrix = zeros((n, m + 1))

                for i, poly in enumerate(polys):
                    for coeff, monom in poly.iter_terms():
                        try:
                            j = list(monom).index(1)
                            matrix[i, j] = coeff
                        except ValueError:
                            matrix[i, m] = -coeff

                soln = solve_linear_system(matrix, *symbols, **flags)
            else:
                soln = solve_poly_system(polys)

            # Use swap_dict to ensure we return the same type as what was
            # passed
            if symbol_swapped:
                if isinstance(soln, dict):
                    res = {}
                    for k in soln.keys():
                        res.update({swap_back_dict[k]: soln[k]})
                    return res
                else:
                    return soln
            else:
                return soln

def solve_linear_system(system, *symbols, **flags):
    """Solve system of N linear equations with M variables, which means
       both Cramer and over defined systems are supported. The possible
       number of solutions is zero, one or infinite. Respectively this
       procedure will return None or dictionary with solutions. In the
       case of over defined system all arbitrary parameters are skipped.
       This may cause situation in with empty dictionary is returned.
       In this case it means all symbols can be assigned arbitrary values.

       Input to this functions is a Nx(M+1) matrix, which means it has
       to be in augmented form. If you are unhappy with such setting
       use 'solve' method instead, where you can input equations
       explicitly. And don't worry about the matrix, this function
       is persistent and will make a local copy of it.

       The algorithm used here is fraction free Gaussian elimination,
       which results, after elimination, in upper-triangular matrix.
       Then solutions are found using back-substitution. This approach
       is more efficient and compact than the Gauss-Jordan method.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       Solve the following system:

              x + 4 y ==  2
           -2 x +   y == 14

       >>> system = Matrix(( (1, 4, 2), (-2, 1, 14)))
       >>> solve_linear_system(system, x, y)
       {x: -6, y: 2}

    """
    matrix = system[:,:]
    syms = list(symbols)

    i, m = 0, matrix.cols-1  # don't count augmentation

    while i < matrix.rows:
        if i == m:
            # an overdetermined system
            if any(matrix[i:,m]):
                return None   # no solutions
            else:
                # remove trailing rows
                matrix = matrix[:i,:]
                break

        if not matrix[i, i]:
            # there is no pivot in current column
            # so try to find one in other columns
            for k in xrange(i+1, m):
                if matrix[i, k]:
                    break
            else:
                if matrix[i, m]:
                    return None   # no solutions
                else:
                    # zero row or was a linear combination of
                    # other rows so now we can safely skip it
                    matrix.row_del(i)
                    continue

            # we want to change the order of colums so
            # the order of variables must also change
            syms[i], syms[k] = syms[k], syms[i]
            matrix.col_swap(i, k)

        pivot_inv = S.One / matrix [i, i]

        # divide all elements in the current row by the pivot
        matrix.row(i, lambda x, _: x * pivot_inv)

        for k in xrange(i+1, matrix.rows):
            if matrix[k, i]:
                coeff = matrix[k, i]

                # subtract from the current row the row containing
                # pivot and multiplied by extracted coefficient
                matrix.row(k, lambda x, j: simplify(x - matrix[i, j]*coeff))

        i += 1

    # if there weren't any problems, augmented matrix is now
    # in row-echelon form so we can check how many solutions
    # there are and extract them using back substitution

    simplified = flags.get('simplified', True)

    if len(syms) == matrix.rows:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k+1, m):
                content -= matrix[k, j]*solutions[syms[j]]

            if simplified:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    elif len(syms) > matrix.rows:
        # this system will have infinite number of solutions
        # dependent on exactly len(syms) - i parameters
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k+1, i):
                content -= matrix[k, j]*solutions[syms[j]]

            # run back-substitution for parameters
            for j in xrange(i, m):
                content -= matrix[k, j]*syms[j]

            if simplified:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    else:
        return None   # no solutions

def solve_undetermined_coeffs(equ, coeffs, sym, **flags):
    """Solve equation of a type p(x; a_1, ..., a_k) == q(x) where both
       p, q are univariate polynomials and f depends on k parameters.
       The result of this functions is a dictionary with symbolic
       values of those parameters with respect to coefficients in q.

       This functions accepts both Equations class instances and ordinary
       SymPy expressions. Specification of parameters and variable is
       obligatory for efficiency and simplicity reason.

       >>> from sympy import *
       >>> a, b, c, x = symbols('a', 'b', 'c', 'x')

       >>> solve_undetermined_coeffs(Eq(2*a*x + a+b, x), [a, b], x)
       {a: 1/2, b: -1/2}

       >>> solve_undetermined_coeffs(Eq(a*c*x + a+b, x), [a, b], x)
       {a: 1/c, b: -1/c}

    """
    if isinstance(equ, Equality):
        # got equation, so move all the
        # terms to the left hand side
        equ = equ.lhs - equ.rhs

    system = collect(equ.expand(), sym, evaluate=False).values()

    if not any([ equ.has(sym) for equ in system ]):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian elimination algorithm
        return solve(system, *coeffs, **flags)
    else:
        return None # no solutions

def solve_linear_system_LU(matrix, syms):
    """ LU function works for invertible only """
    assert matrix.rows == matrix.cols-1
    A = matrix[:matrix.rows,:matrix.rows]
    b = matrix[:,matrix.cols-1:]
    soln = A.LUsolve(b)
    solutions = {}
    for i in range(soln.rows):
        solutions[syms[i]] = soln[i,0]
    return solutions

def dsolve(eq, funcs):
    """
    Solves any (supported) kind of differential equation.

    Usage

        dsolve(f, y(x)) -> Solve a differential equation f for the function y


    Details

        @param f: ordinary differential equation (either just the left hand
            side, or the Equality class)

        @param y: indeterminate function of one variable

        - You can declare the derivative of an unknown function this way:
        >>> from sympy import *
        >>> x = Symbol('x') # x is the independent variable
        >>> f = Function("f")(x) # f is a function of x
        >>> f_ = Derivative(f, x) # f_ will be the derivative of \
        f with respect to x

        - This function just parses the equation "eq" and determines the type of
        differential equation by its order, then it determines all the
        coefficients and then calls the particular solver, which just accepts
        the coefficients.
        - "eq" can be either an Equality, or just the left hand side (in which
          case the right hand side is assumed to be 0)
        - See test_ode.py for many tests, which serve also as a set of examples
          how to use dsolve
        - dsolve always returns an equality class.  If possible, it solves the
          solution explicitly for the function being solved for. Otherwise, it
          returns an implicit solution.
        - Arbitrary constants are symbols named C1, C2, and so on.

    Examples

        >>> from sympy import *
        >>> x = Symbol('x')

        >>> f = Function('f')
        >>> dsolve(Derivative(f(x),x,x)+9*f(x), f(x))
        f(x) == C1*sin(3*x) + C2*cos(3*x)
        >>> dsolve(Eq(Derivative(f(x),x,x)+9*f(x)+1, 1), f(x))
        f(x) == C1*sin(3*x) + C2*cos(3*x)

    """

    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return dsolve(eq.lhs-eq.rhs, funcs)
        eq = eq.lhs

    #currently only solve for one function
    if isinstance(funcs, Basic) or len(funcs) == 1:
        if isinstance(funcs, (list, tuple)): # normalize args
            f = funcs[0]
        else:
            f = funcs

        if len(f.args) != 1:
            raise NotImplementedError("Only functions of one variable are supported")
        x = f.args[0]
        f = f.func

        # Collect diff(f(x),x) terms so that match will work correctly
        eq = collect(eq, f(x).diff(x))
        #We first get the order of the equation, so that we can choose the
        #corresponding methods. Currently, only first and second
        #order odes can be handled.
        order = deriv_degree(eq, f(x))

        if order > 1:
            return constantsimp(solve_ODE_higher_order(eq, f(x), order), x, 2*order)
        if  order > 2:
           raise NotImplementedError("dsolve: Cannot solve " + str(eq))
        elif order == 2:
            return solve_ODE_second_order(eq, f(x))
        elif order == 1:
           return constantsimp(solve_ODE_first_order(eq, f(x)), x, 1)
        else:
            raise NotImplementedError("Not a differential equation.")

def classify_ode(expr, func):
    """
    Returns a tuple of possible dsolve classifications.

    The first item in the tuple is the classification that dsolve uses to solve
    the ode by default.  To make dsolve use a different classification, use
    dsolve(ode, func, hint=<classification>).  To make dsolve apply all relevant
    classification hints, use dsolve(ode, func, hint="all").  This will return
    a dictionary of each method with its corresponding solution.  To have dsolve
    try all methods and return the simplest one, use dsolve(ode, func,
    hint="best").   This takes into account whether the solution is solvable in
    the function, whether it contains any Integral classes, and which one is the
    shortest in size.

    Some classifications are duplicated with the word _alt at the end.  This
    means that the dsolve method for that classification can produce more than
    one solution.  Use the _alt classification to obtain the non-default
    solution.  For example, first order differential equations with homogeneous
    coefficients will return both "1st_homogeneous_coeff" and
    "1st_homogeneous_coeff_alt".

    Because all solutions should be mathematically equivalent, some dsolve hints
    may return the exact same result for an ode.  Often, though, two different
    hints will return the same solution formatted differently.  The two should
    be equivalent.
    """
    # A list of hints in the order that they should be applied.  That means
    # that, in general, hints earlier in the list should produce simpler results
    # than those later for odes that fit both.  This is just based on my own
    # empirical observations, so if you find that *in general*, a hint later in
    # the list is better than one before it, fell free to modify the list.  Note
    # however that you can easily override in dsolve (see the docstring).
    hints = ("separable", "exact", "1st_linear", "Bernoulli",\
    "1st_homogeneous_coeff", "1st_homogeneous_coeff_alt", "nth_homogeneous")

    # This will have to wait for nth order homogeneous because I will need to
    # clean up solve_ODE_second_order and solve_ODE_1 first.

@vectorize(0)
def constantsimp(expr, independentsymbol, endnumber, startnumber=1,
    symbolname='C'):
    """
    Simplifies an expression with arbitrary constants in it.

    This function is written specifically to work with dsolve(), and is not
    indented for general use.

    This is done by "absorbing" the arbitrary constants in to other arbitrary
    constants, numbers, and symbols for which they are not independent of.

    The symbols must all have the same name with numbers after it, for example,
    C1, C2, C3.  The symbolname here would be 'C', the startnumber would be 1,
    and the end number would be 3.  If the arbitrary constants are independent
    of the variable x, then the independentsymbol would be x.

    Because terms are "absorbed" into arbitrary constants and because constants
    are renumbered after simplifying, the arbitrary constants in expr are not
    necessarily equal to the ones of the same name in the returned result.

    If two or more arbitrary constants are added, multiplied, or raised to the
    power of each other, they are first absorbed together into a single
    arbitrary constant.  Then the new constant is combined into other terms
    if necessary.

    Absorption is done naively.  constantsimp() does not attempt to expand
    or simplify the expression first to obtain better absorption.

    Constants are renumbered after simplification so that they are sequential,
    such as C1, C2, C3, and so on.

    Example:
    >>> from sympy import *
    >>> C1, C2, C3, x, y = symbols('C1 C2 C3 x y')
    >>> constantsimp(2*C1*x, x, 3)
    C1*x
    >>> constantsimp(C1 + 2 + x + y, x, 3)
    C1 + x
    >>> constantsimp(C1*C2 + 2 + x + y + C3*x, x, 3)
    C1 + x + C2*x
    """
    # We need to have an internal recursive function so that newstartnumber
    # maintains its values throughout recursive calls

    global newstartnumber
    newstartnumber = 1

    def _constantsimp(expr, independentsymbol, endnumber, startnumber=1,
    symbolname='C'):
        """
        The function works recursively.  The idea is that, for Mul, Add, Pow, and
        Function, if the class has a constant in it, then we can simplify it,
        which we do by recursing down and simplifying up.  Otherwise, we can skip
        that part of the expression.
        """
        constantsymbols = [Symbol(symbolname+"%d" % t) for t in range(startnumber,
        endnumber + 1)]
        x = independentsymbol

        if isinstance(expr, Equality):
            return Equality(_constantsimp(expr.lhs, x, endnumber, startnumber,
                symbolname), _constantsimp(expr.rhs, x, endnumber, startnumber,
                symbolname))

        if type(expr) not in (Mul, Add, Pow) and not expr.is_Function:
            # We don't know how to handle other classes
            # This also serves as the base case for the recursion
            return expr
        elif not any(t in expr for t in constantsymbols):
            return expr
        else:
            newargs = []
            hasconst = False
            isPowExp = False
            reeval = False
            for i in expr.args:
                if i not in constantsymbols:
                    newargs.append(i)
                else:
                    newconst = i
                    hasconst = True
                    if expr.is_Pow and i == expr.exp:
                        isPowExp = True

            for i in range(len(newargs)):
                isimp = _constantsimp(newargs[i], x, endnumber, startnumber,
                symbolname)
                if isimp in constantsymbols:
                    reeval = True
                    hasconst = True
                    newconst = isimp
                    if expr.is_Pow and i == 1:
                        isPowExp = True
                newargs[i] = isimp
            if hasconst:
                newargs = [i for i in newargs if i.has(x)]
                if isPowExp:
                    newargs = newargs + [newconst] # Order matters in this case
                else:
                    newargs = [newconst] + newargs
            if expr.is_Pow and len(newargs) == 1:
                newargs.append(S.One)
            if expr.is_Function:
                if (len(newargs) == 0 or hasconst and len(newargs) == 1):
                    return newconst
                else:
                    newfuncargs = [_constantsimp(t, x, endnumber, startnumber,
                    symbolname) for t in expr.args]
                    return expr.new(*newfuncargs)
            else:
                newexpr = expr.new(*newargs)
                if reeval:
                    return _constantsimp(newexpr, x, endnumber, startnumber,
                    symbolname)
                else:
                    return newexpr

    def _renumber(expr, symbolname, startnumber, endnumber):
        """
        Renumber arbitrary constants in expr.

        This is a simple function that goes through and renumbers any Symbol
        with a name in the form symbolname + num where num is in the range
        from startnumber to endnumber.

        Symbols are renumbered in the order that they are encountered via
        a depth first search through args, so they should be numbered roughly
        in the order that they appear in the final, printed expression.

        The structure of the function is very similar to _constantsimp().
        """
        constantsymbols = [Symbol(symbolname+"%d" % t) for t in range(startnumber,
        endnumber + 1)]
        global newstartnumber

        if isinstance(expr, Equality):
            return Equality(_renumber(expr.lhs, symbolname, startnumber, endnumber),
            _renumber(expr.rhs, symbolname, startnumber, endnumber))

        if type(expr) not in (Mul, Add, Pow) and not expr.is_Function and\
        not any(t in expr for t in constantsymbols):
            # Base case, as above.  We better hope there aren't constants inside
            # of some other class, because they won't be simplified.
            return expr
        elif expr in constantsymbols:
            # Renumbering happens here
            newconst = Symbol(symbolname + str(newstartnumber))
            newstartnumber += 1
            return newconst
        else:
            sortedargs = list(expr.args)
            sortedargs.sort(Basic._compare_pretty)
            if expr.is_Function or expr.is_Pow:
                return expr.new(*map(lambda x: _renumber(x, symbolname, \
                startnumber, endnumber), expr.args))
            else:
                return expr.new(*map(lambda x: _renumber(x, symbolname, \
                startnumber, endnumber), sortedargs))


    simpexpr = _constantsimp(expr, independentsymbol, endnumber, startnumber,
    symbolname)

    return _renumber(simpexpr, symbolname, startnumber, endnumber)



def deriv_degree(expr, func):
    """ get the order of a given ode, the function is implemented
    recursively """
    a = Wild('a', exclude=[func])

    order = 0
    if isinstance(expr, Derivative):
        order = len(expr.symbols)
    else:
        for arg in expr.args:
            if isinstance(arg, Derivative):
                order = max(order, len(arg.symbols))
            elif expr.match(a):
                order = 0
            else :
                for arg1 in arg.args:
                    order = max(order, deriv_degree(arg1, func))

    return order

def solve_ODE_first_order(eq, f):
    """
    Solves many kinds of first order odes.
    Different methods are used depending on the form of the given equation.
    Now the linear, Bernoulli, exact, and first order homogeneous cases are
    implemented.
    """
    from sympy.integrals.integrals import integrate
    x = f.args[0]
    f = f.func
    C1 = Symbol('C1')

    # Linear case: a(x)*y'+b(x)*y+c(x) == 0
    a = Wild('a', exclude=[f(x)])
    b = Wild('b', exclude=[f(x)])
    c = Wild('c', exclude=[f(x)])

    r = eq.match(a*diff(f(x),x) + b*f(x) + c)
    if r:
        t = exp(integrate(r[b]/r[a], x))
        tt = integrate(t*(-r[c]/r[a]), x)
        return Equality(f(x),(tt + C1)/t)

    # Bernoulli case: a(x)*y'+b(x)*y+c(x)*y**n == 0
    n = Wild('n', exclude=[f(x)])

    r = eq.match(a*diff(f(x),x) + b*f(x) + c*f(x)**n)

    if r:
        if r[n] != 1:
            t = C.exp((1-r[n])*integrate(r[b]/r[a],x))
            tt = (r[n]-1)*integrate(t*r[c]/r[a],x)
            return Equality(f(x),((tt + C1)/t)**(1/(1-r[n])))
        #if r[n] == 1:
         #   return Equality(f(x),C1*exp(integrate(-(r[b]+r[c]), x)))

    a = Wild('a', exclude=[f(x).diff(x)])
    b = Wild('b', exclude=[f(x).diff(x)])
    r = eq.match(a+b*diff(f(x),x))
    # This match is used for several cases below.
    if r:
        y = Symbol('y', dummy=True)
        r[a] = r[a].subs(f(x),y)
        r[b] = r[b].subs(f(x),y)

        # Separable Case: y' == P(y)*Q(x)
        r[a] = separatevars(r[a])
        r[b] = separatevars(r[b])
        # m1[coeff]*m1[x]*m1[y] + m2[coeff]*m2[x]*m2[y]*y'
        m1 = separatevars(r[a], dict=True, symbols=(x, y))
        m2 = separatevars(r[b], dict=True, symbols=(x, y))

        if m1 and m2:
            return Equality(integrate(m2['coeff']*m2[y]/m1[y], y).subs(y, f(x)),\
            integrate(-m1['coeff']*m1[x]/m2[x], x)+C1)
        # Exact Differential Equation: P(x,y)+Q(x,y)*y'=0 where dP/dy == dQ/dx
        if simplify(r[a].diff(y)) == simplify(r[b].diff(x)) and r[a]!=0:
            x0 = Symbol('x0', dummy=True)
            y0 = Symbol('y0', dummy=True)
            tmpsol = integrate(r[b].subs(x,x0),(y,y0,y))+integrate(r[a],(x,x0,x))
            sol = 0
            assert tmpsol.is_Add
            for i in tmpsol.args:
                if x0 not in i and y0 not in i:
                    sol += i
            assert sol != 0
            sol = Equality(sol,C1)

            try:
                # See if the equation can be solved explicitly for f
                # This part of the code will change when solve returns RootOf.
                sol1 = solve(sol,y)
            except NotImplementedError:
                return sol.subs(y,f(x))
            else:
                if len(sol1) !=1:
                    return sol.subs(y,f(x))
                else:
                    return Equality(f(x),sol1[0].subs(y,f(x)))

        # First order equation with homogeneous coefficients:
        # dy/dx == F(y/x) or dy/dx == F(x/y)
        ordera = homogeneous_order(r[a], x, y)
        orderb = homogeneous_order(r[b], x, y)
        if ordera == orderb and ordera != None:
            # There are two substitutions that solve the equation, u=x/y and u=y/x
            # They produce different integrals, so try them both and see which
            # one is easier.
            u1 = Symbol('u1', dummy=True) # u1 == y/x
            u2 = Symbol('u2', dummy=True) # u2 == x/y
            _a = Symbol('_a', dummy=True)
            #print ((-r[b]/(r[a]+u1*r[b])).subs({x:1, y:u1}),
            #print (-r[a]/(r[b]+u2*r[a])).subs({x:u2, y:1}))
            int1 = integrate((-r[b]/(r[a]+u1*r[b])).subs({x:1, y:u1}), u1)
            int2 = integrate((-r[a]/(r[b]+u2*r[a])).subs({x:u2, y:1}), u2)
            # Substitute back in for u1 and u2.
            if int1.has(C.Integral):
                int1 = C.Integral(int1.args[0],(u1,_a,f(x)/x))
            else:
                int1 = int1.subs(u1,f(x)/x)
            if int2.has(C.Integral):
                int2 = C.Integral(int2.args[0],(u2,_a,x/f(x)))
            else:
                int2 = int2.subs(u2,x/f(x))
            sol1 = logcombine(Equality(log(x), int1 + log(C1)), assume_pos_real=True)
            sol2 = logcombine(Equality(log(f(x)), int2 + log(C1)), assume_pos_real=True)
            if sol1.lhs.is_Function and sol1.lhs.func == log and sol1.rhs == 0:
                sol1 = Equality(sol1.lhs.args[0]*C1,C1)
            if sol2.lhs.is_Function and sol2.lhs.func == log and sol2.rhs == 0:
                sol2 = Equality(sol2.lhs.args[0]*C1,C1)

            # There are two solutions.  We need to determine which one to use
            # First, if they are the same, don't bother testing which one to use
            if sol1 == sol2:
                # But still try to solve for f
                try:
                    sol1s = map((lambda t: t.subs(y, f(x))),\
                    solve(sol1.lhs.subs(f(x),y)-sol1.rhs.subs(f(x),y), y))
                    if sol1s == []:
                        raise NotImplementedError
                except NotImplementedError:
                    return sol1
                else:
                    sol1sr = map((lambda t: Equality(f(x), t.subs({u1:f(x)/x,\
                    y:f(x)}))), sol1s)
                    if len(sol1sr) == 1:
                        return logcombine(sol1sr[0], assume_pos_real=True)
                    else:
                        return [logcombine(t, assume_pos_real=True) for t in sol1sr]
            # Second, try to return an evaluated integral:
            if sol1.has(C.Integral):
                return sol2
            if sol2.has(C.Integral):
                return sol1
            # Next, try to return an explicit solution.  This code will change
            # when RootOf is implemented in Solve.
            try:
                sol1s = map((lambda t: t.subs(y, f(x))),\
                solve(sol1.lhs.subs(f(x),y)-sol1.rhs.subs(f(x),y), y))
                if sol1s == []:
                    raise NotImplementedError
            except NotImplementedError:
                pass
            else:
                sol1sr = map((lambda t: Equality(f(x), t.subs({u1:f(x)/x,\
                y:f(x)}))), sol1s)
                if len(sol1sr) == 1:
                    return logcombine(sol1sr[0], assume_pos_real=True)
                else:
                    return [logcombine(t, assume_pos_real=True) for t in sol1sr]
            try:
                sol2s = map((lambda t: t.subs(y, f(x))),\
                solve(sol2.lhs.subs(f(x),y)-sol2.rhs.subs(f(x),y), y))
                if sol2s == []:
                    raise NotImplementedError
            except NotImplementedError:
                pass
            else:
                sol2sr = map((lambda t: Equality(f(x), t.subs({u2:x/f(x),\
                y:f(x)}))), sol2s)
                if len(sol2sr) == 1:
                    return logcombine(sol2sr[0], assume_pos_real=True)
                else:
                    return [logcombine(t, assume_pos_real=True) for t in sol2srs]

            # Finally, try to return the shortest expression, naively computed
            # based on the length of the string version of the expression.  This
            # may favor combined fractions because they will not have duplicate
            # denominators, and may slightly favor expressions with fewer
            # additions and subtractions, as those are separated by spaces by
            # the printer.
            return min(sol1, sol2, key=(lambda x: len(str(x))))

    # Other cases of first order odes will be implemented here

    raise NotImplementedError("solve_ODE_first_order: Cannot solve " + str(eq))

def solve_ODE_higher_order(eq, f, order):
    from sympy.integrals.integrals import integrate
    x = f.args[0]
    f = f.func
    b = Wild('b', exclude=[f(x)])
    j = 0
    s = S(0)
    wilds = []
    constants = numbered_symbols(prefix='C', function=Symbol, start=1)
    for i in numbered_symbols(prefix='a', function=Wild, exclude=[f(x)]):
        if j == order+1:
            break
        wilds.append(i)
        s += i*f(x).diff(x,j)
        j += 1
    s += b

    r = eq.match(s)
    if r:
        # The ODE is homogeneous
        if all([not r[i].has(x) for i in wilds]):
            # First, set up characteristic equation.
            m = Symbol('m', dummy=True)
            chareq = S(0)
            for i in r:
                if i == b:
                    pass
                else:
                    chareq += r[i]*m**S(i.name[1:])
            chareqroots = RootsOf(chareq, m)
            charroots_exact = list(chareqroots.exact_roots())
            charroots_formal = list(chareqroots.formal_roots())
            if charroots_formal and discriminant(chareq, m) == 0:
                # If Poly cannot find the roots explicitly, we can only return
                # an expression in terms of RootOf's if we know the roots
                # are not repeated.  We use the fact that a polynomial has
                # repeated roots iff its discriminant == 0.

                # TODO: cancel out roots from charroots_exact, then check
                # the discriminant of chareq.
                raise NotImplementedError("Cannot find all of the roots of " + \
                "characteristic equation " + str(chareq) + ", which has " + \
                "repeated roots.")
            # Create a dict root: multiplicity or charroots
            charroots = {}
            for i in charroots_exact + charroots_formal:
                if i in charroots:
                    charroots[i] += 1
                else:
                    charroots[i] = 1
            gsol = S(0)
            psol = S(0)
            # We need keep track of terms so we can run collect() at the end.
            # This is necessary for constantsimp to work properly.
            collectterms = []
            for root, multiplicity in charroots.items():
                for i in range(multiplicity):
                    if isinstance(root, RootOf):
                        # re and im do not work with RootOf, so the work around is
                        # to put solution in non (complex) expanded form.
                        # See issue 1563.
                        gsol += exp(root*x)*constants.next()
                        assert multiplicity == 1
                    else:
                        reroot = re(root)
                        imroot = im(root)
                        gsol += x**i*exp(reroot*x)*(constants.next()*sin(abs(imroot)*x) \
                        + constants.next()*cos(imroot*x))
                        collectterms = [(i, reroot, imroot)] + collectterms
            gsol = expand_mul(gsol, deep=False)
            for i, reroot, imroot in collectterms:
                gsol = collect(gsol, x**i*exp(reroot*x)*sin(abs(imroot)*x))
                gsol = collect(gsol, x**i*exp(reroot*x)*cos(imroot*x))
            for i, reroot, imroot in collectterms:
                gsol = collect(gsol, x**i*exp(reroot*x))
            if r[b] != 0:
                # Variation of Paramters
                gensols = []
                # Keep track of when to use sin or cos for nonzero imroot
                trigdict = {}
                if len(collectterms) != order:
                    raise NotImplementedError("Cannot find " + str(order) + \
                    " solutions to homogeneous equation to apply variation " + \
                    "of parameters to " + str(eq))
                for i, reroot, imroot in collectterms:
                    if imroot == 0:
                        gensols.append(x**i*exp(reroot*x))
                    else:
                        if x**i*exp(reroot*x)*sin(abs(imroot)*x) in gensols:
                            gensols.append(x**i*exp(reroot)*cos(imroot*x))
                        else:
                            gensols.append(x**i*exp(reroot*x)*sin(abs(imroot)*x))
                wr = wronskian(gensols, x)
                wr = trigsimp(wr) # to reduce sin(x)**2 + cos(x)**2 to 1
                print 'wr', wr
                if not wr:
                    raise NotImplementedError("Cannot find " + str(order) + \
                    " solutions to homogeneous equation to apply variation " + \
                    "of parameters to " + str(eq) + " (Wronskian == 0)")
                negoneterm = (-1)**(order)
                for i in gensols:
                    psol += negoneterm*C.Integral(wronskian(filter(lambda x: x != i, \
                    gensols), x)*r[b]/wr, x)*i/r[wilds[-1]]
                    negoneterm *= -1
            psol = simplify(psol)
            psol = trigsimp(psol, deep=True)
            from sympy import pprint
            pprint(psol)
            return Equality(f(x), gsol + psol)



    # special equations, that we know how to solve
    # TODO: refactor into substitution u = f' (divide out exp(-f(x))
    a = Wild('a')
    t = x*exp(f(x))
    tt = a*t.diff(x, x)/t
    r = eq.match(tt.expand())
    if r:
        return Equality(f(x),log(constants.next()+constants.next()/x))

    t = x*exp(-f(x))
    tt = a*t.diff(x, x)/t
    r = eq.match(tt.expand())
    if r:
        #check, that we've rewritten the equation correctly:
        #assert ( r[a]*t.diff(x,2)/t ) == eq.subs(f, t)
        return Equality(f(x),-log(constants.next()+constants.next()/x))

    neq = eq*exp(f(x))/exp(-f(x))
    r = neq.match(tt.expand())
    if r:
        #check, that we've rewritten the equation correctly:
        #assert ( t.diff(x,2)*r[a]/t ).expand() == eq
        return Equality(f(x),-log(constants.next()+constants.next()/x))


    raise NotImplementedError("solve_ODE_higher_order: Cannot solve " + str(eq)) # Yet!



def solve_ODE_second_order(eq, f):
    """
    solves many kinds of second order odes, different methods are used
    depending on the form of the given equation. So far the constants
    coefficients case and a special case are implemented.
    """
    x = f.args[0]
    f = f.func
    C1 = Symbol('C1')
    C2 = Symbol('C2')

    #constant coefficients case: af''(x)+bf'(x)+cf(x)=0
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])

    r = eq.match(a*f(x).diff(x,x) + c*f(x))
    if r:
        return Equality(f(x),C1*C.sin(sqrt(r[c]/r[a])*x)+C2*C.cos(sqrt(r[c]/r[a])*x))

    r = eq.match(a*f(x).diff(x,x) + b*diff(f(x),x) + c*f(x))
    if r:
        r1 = solve(r[a]*x**2 + r[b]*x + r[c], x)
        if r1[0].is_real:
            if len(r1) == 1:
                return Equality(f(x),(C1 + C2*x)*exp(r1[0]*x))
            else:
                return Equality(f(x),C1*exp(r1[0]*x) + C2*exp(r1[1]*x))
        else:
            r2 = abs((r1[0] - r1[1])/(2*S.ImaginaryUnit))
            return Equality(f(x),(C2*C.cos(r2*x) + C1*C.sin(r2*x))*exp((r1[0] + r1[1])*x/2))

    #other cases of the second order odes will be implemented here



def homogeneous_order(eq, *symbols):
    """
    Determines if a function is homogeneous and if so of what order.
    A function f(x,y,...) is homogeneous of order n if
    f(t*x,t*y,t*...) == t**n*f(x,y,...).  It is implemented recursively.

    Functions can be symbols, but every argument of the function must also be
    a symbol, and the arguments of the function that appear in the expression
    must match those given in the list of symbols.  If a declared function
    appears with different arguments than given in the list of symbols, None is
    returned.
    Example:
    >>> from sympy import *
    >>> from sympy.solvers.solvers import homogeneous_order
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> f = Function('f')
    >>> homogeneous_order(f(x), f(x)) == None
    True
    >>> homogeneous_order(f(x,y), f(y, x), x, y) == None
    True
    >>> homogeneous_order(f(x), f(x), x)
    1

    Returns the order n if g is homogeneous and None if it is not homogeneous.
    Examples:
    >>> homogeneous_order(x**2*f(x)/sqrt(x**2+f(x)**2), x, f(x))
    2
    >>> homogeneous_order(x**2+f(x), x, f(x)) == None
    True
    """
    if eq.has(log):
        eq = logcombine(eq, assume_pos_real=True)
    # This runs as a separate function call so that logcombine doesn't endlessly
    # put back together what homogeneous_order is trying to take apart.
    return _homogeneous_order(eq, *symbols)

def _homogeneous_order(eq, *symbols):

    if not symbols:
        raise ValueError, "homogeneous_order: no symbols were given."

    n = set()

    # Replace all functions with dummy variables

    if any(getattr(i, 'is_Function') for i in symbols):
        for i in symbols:
            if i.is_Function:
                if not all(map((lambda i: i in symbols), i.args)):
                    return None
                elif i not in symbols:
                    pass
                else:
                    dummyvar = numbered_symbols(prefix='d', dummy=True).next()
                    eq = eq.subs(i, dummyvar)
                    symbols = list(symbols)
                    symbols.remove(i)
                    symbols.append(dummyvar)
                    symbols = tuple(symbols)

    # The following are not supported
    if eq.is_Order or eq.is_Derivative:
        return None

    # These are all constants
    if type(eq) in (int, float) or eq.is_Number or eq.is_Integer or \
    eq.is_Rational or eq.is_NumberSymbol or eq.is_Real:
        return sympify(0)

    # Break the equation into additive parts
    if eq.is_Add:
        s = set()
        for i in eq.args:
            s.add(_homogeneous_order(i, *symbols))
        if len(s) != 1:
            return None
        else:
            n = s

    if eq.is_Pow:
        if not eq.args[1].is_Number:
            return None
        o = _homogeneous_order(eq.args[0], *symbols)
        if o == None:
            return None
        else:
            n.add(sympify(o*eq.args[1]))

    t = Symbol('t', dummy=True, positive=True) # It is sufficient that t > 0
    r = Wild('r', exclude=[t])
    a = Wild('a', exclude=[t])
    eqs = eq.subs(dict(zip(symbols,(t*i for i in symbols))))

    if eqs.is_Mul:
        if t not in eqs:
            n.add(sympify(0))
        else:
            m = eqs.match(r*t**a)
            if m:
                n.add(sympify(m[a]))
            else:
                s = 0
                for i in eq.args:
                    o = _homogeneous_order(i, *symbols)
                    if o == None:
                        return None
                    else:
                        s += o
                n.add(sympify(s))

    if eq.is_Function:
        if eq.func == log:
            # The only possibility to pull a t out of a function is a power in
            # a logarithm.  This is very likely due to calling of logcombine.
            if eq.args[0].is_Pow:
                return _homogeneous_order(eq.args[0].args[1]*log(eq.args[0].args[0]), *symbols)
            elif eq.args[0].is_Mul and all(i.is_Pow for i in iter(eq.args[0].args)):
                arg = 1
                pows = set()
                for i in eq.args[0].args:
                    if i.args[1].args[0] == -1:
                        arg *= 1/i.args[0]
                        pows.add(sympify(-1*i.args[1]))
                    else:
                        arg *= i.args[0]
                        pows.add(sympify(i.args[1]))
                if len(pows) != 1:
                    return None
                else:
                    return _homogeneous_order(pows.pop()*log(arg), *symbols)
            else:
                if _homogeneous_order(eq.args[0], *symbols) == 0:
                    return sympify(0)
                else:
                    return None
        else:
            if _homogeneous_order(eq.args[0], *symbols) == 0:
                return sympify(0)
            else:
                return None

    if len(n) != 1 or n == None:
        return None
    else:
        return n.pop()

    return None

x = Symbol('x', dummy=True)
a,b,c,d,e,f,g,h = [Wild(t, exclude=[x]) for t in 'abcdefgh']
patterns = None

def _generate_patterns():
    """
    Generates patterns for transcendental equations.

    This is lazily calculated (called) in the tsolve() function and stored in
    the patterns global variable.
    """

    tmp1 = f ** (h-(c*g/b))
    tmp2 = (-e*tmp1/a)**(1/d)
    global patterns
    patterns = [
        (a*(b*x+c)**d + e   , ((-(e/a))**(1/d)-c)/b),
        (    b+c*exp(d*x+e) , (log(-b/c)-e)/d),
        (a*x+b+c*exp(d*x+e) , -b/a-LambertW(c*d*exp(e-b*d/a)/a)/d),
        (    b+c*f**(d*x+e) , (log(-b/c)-e*log(f))/d/log(f)),
        (a*x+b+c*f**(d*x+e) , -b/a-LambertW(c*d*f**(e-b*d/a)*log(f)/a)/d/log(f)),
        (    b+c*log(d*x+e) , (exp(-b/c)-e)/d),
        (a*x+b+c*log(d*x+e) , -e/d+c/a*LambertW(a/c/d*exp(-b/c+a*e/c/d))),
        (a*(b*x+c)**d + e*f**(g*x+h) , -c/b-d*LambertW(-tmp2*g*log(f)/b/d)/g/log(f))
    ]

def tsolve(eq, sym):
    """
    Solves a transcendental equation with respect to the given
    symbol. Various equations containing mixed linear terms, powers,
    and logarithms, can be solved.

    Only a single solution is returned. This solution is generally
    not unique. In some cases, a complex solution may be returned
    even though a real solution exists.

        >>> from sympy import *
        >>> x = Symbol('x')

        >>> tsolve(3**(2*x+5)-4, x)
        [(-5*log(3) + log(4))/(2*log(3))]

        >>> tsolve(log(x) + 2*x, x)
        [LambertW(2)/2]

    """
    if patterns is None:
        _generate_patterns()
    eq = sympify(eq)
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    sym = sympify(sym)
    eq2 = eq.subs(sym, x)
    # First see if the equation has a linear factor
    # In that case, the other factor can contain x in any way (as long as it
    # is finite), and we have a direct solution to which we add others that
    # may be found for the remaining portion.
    r = Wild('r')
    m = eq2.match((a*x+b)*r)
    if m and m[a]:
        return [(-b/a).subs(m).subs(x, sym)] + solve(m[r], x)
    for p, sol in patterns:
        m = eq2.match(p)
        if m:
            return [sol.subs(m).subs(x, sym)]

    # let's also try to inverse the equation
    lhs = eq
    rhs = S.Zero

    while True:
        indep, dep = lhs.as_independent(sym)

        # dep + indep == rhs
        if lhs.is_Add:
            # this indicates we have done it all
            if indep is S.Zero:
                break

            lhs = dep
            rhs-= indep

        # dep * indep == rhs
        else:
            # this indicates we have done it all
            if indep is S.One:
                break

            lhs = dep
            rhs/= indep

    #                    -1
    # f(x) = g  ->  x = f  (g)
    if lhs.is_Function and lhs.nargs==1 and hasattr(lhs, 'inverse'):
        rhs = lhs.inverse() (rhs)
        lhs = lhs.args[0]

        sol = solve(lhs-rhs, sym)
        return sol

    elif lhs.is_Add:
        # just a simple case - we do variable substitution for first function,
        # and if it removes all functions - let's call solve.
        #      x    -x                   -1
        # UC: e  + e   = y      ->  t + t   = y
        t = Symbol('t', dummy=True)
        terms = lhs.args

        # find first term which is Function
        for f1 in lhs.args:
            if f1.is_Function:
                break
        else:
            assert False, 'tsolve: at least one Function expected at this point'

        # perform the substitution
        lhs_ = lhs.subs(f1, t)

        # if no Functions left, we can proceed with usual solve
        if not (lhs_.is_Function or
                any(term.is_Function for term in lhs_.args)):
            cv_sols = solve(lhs_ - rhs, t)
            for sol in cv_sols:
                if sol.has(sym):
                    raise NotImplementedError("Unable to solve the equation")
            cv_inv = solve( t - f1, sym )[0]
            sols = list()
            for sol in cv_sols:
                sols.append(cv_inv.subs(t, sol))
            return sols


    raise NotImplementedError("Unable to solve the equation.")

def msolve(*args, **kwargs):
    """
    Compatibility wrapper pointing to nsolve().

    msolve() has been renamed to nsolve(), please use nsolve() directly."""
    warn('msolve() is has been renamed, please use nsolve() instead',
         DeprecationWarning)
    args[0], args[1] = args[1], args[0]
    return nsolve(*args, **kwargs)

# TODO: option for calculating J numerically
def nsolve(*args, **kwargs):
    """
    Solve a nonlinear equation system numerically.

    nsolve(f, [args,] x0, modules=['mpmath'], **kwargs)

    f is a vector function of symbolic expressions representing the system.
    args are the variables. If there is only one variable, this argument can be
    omitted.
    x0 is a starting vector close to a solution.

    Use the modules keyword to specify which modules should be used to evaluate
    the function and the Jacobian matrix. Make sure to use a module that
    supports matrices. For more information on the syntax, please see the
    docstring of lambdify.

    Overdetermined systems are supported.

    >>> from sympy import Symbol, nsolve
    >>> import sympy
    >>> sympy.mpmath.mp.dps = 15
    >>> x1 = Symbol('x1')
    >>> x2 = Symbol('x2')
    >>> f1 = 3 * x1**2 - 2 * x2**2 - 1
    >>> f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
    >>> print nsolve((f1, f2), (x1, x2), (-1, 1))
    [-1.19287309935246]
    [ 1.27844411169911]

    For one-dimensional functions the syntax is simplified:

    >>> from sympy import sin
    >>> nsolve(sin(x), x, 2)
    3.14159265358979
    >>> nsolve(sin(x), 2)
    3.14159265358979

    mpmath.findroot is used, you can find there more extensive documentation,
    especially concerning keyword parameters and available solvers.
    """
    # interpret arguments
    if len(args) == 3:
        f = args[0]
        fargs = args[1]
        x0 = args[2]
    elif len(args) == 2:
        f = args[0]
        fargs = None
        x0 = args[1]
    elif len(args) < 2:
        raise TypeError('nsolve expected at least 2 arguments, got %i'
                        % len(args))
    else:
        raise TypeError('nsolve expected at most 3 arguments, got %i'
                        % len(args))
    modules = kwargs.get('modules', ['mpmath'])
    if isinstance(f,  (list,  tuple)):
        f = Matrix(f).T
    if not isinstance(f, Matrix):
        # assume it's a sympy expression
        if isinstance(f, Equality):
            f = f.lhs - f.rhs
        f = f.evalf()
        atoms = set(s for s in f.atoms() if isinstance(s, Symbol))
        if fargs is None:
            fargs = atoms.copy().pop()
        if not (len(atoms) == 1 and (fargs in atoms or fargs[0] in atoms)):
            raise ValueError('expected a one-dimensional and numerical function')
        f = lambdify(fargs, f, modules)
        return findroot(f, x0, **kwargs)
    if len(fargs) > f.cols:
        raise NotImplementedError('need at least as many equations as variables')
    verbose = kwargs.get('verbose', False)
    if verbose:
        print 'f(x):'
        print f
    # derive Jacobian
    J = f.jacobian(fargs)
    if verbose:
        print 'J(x):'
        print J
    # create functions
    f = lambdify(fargs, f.T, modules)
    J = lambdify(fargs, J, modules)
    # solve the system numerically
    x = findroot(f, x0, J=J, **kwargs)
    return x
