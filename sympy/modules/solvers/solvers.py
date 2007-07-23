
""" This module contain solvers for all kinds of equations:

    - algebraic, use solve()

    - recurrence, use rsolve() (not implemented)

    - differential, use dsolve() (not implemented)

"""

from sympy import *

from sympy.modules.utilities import any
from sympy.modules.matrices import zeronm
from sympy.modules.polynomials import roots
from sympy.modules.simplify import simplify, collect

### NOTE: set simplified=True when 'simplify' module will be merged !!!

def solve(eq, syms, simplified=False):
    """Solves univariate polynomial equations and linear systems with
       arbitrary symbolic coefficients. This function is just a wrapper
       which makes analysis of its arguments and executes more specific
       functions like 'roots' or 'solve_linear_system' etc.

       On input you have to specify equation or a set of equations
       (in this case via a list) using '==' pretty syntax or via
       ordinary expressions, and a list of variables.

       On output you will get a list of solutions in univariate case
       or a dictionary with variables as keys and solutions as values
       in the other case. If there were variables with can be assigned
       with arbitrary value, then they will be avoided in the output.

       Optionaly it is possible to have the solutions preprocessed
       using simplification routines if 'simplified' flag is set.

       To solve recurrence relations or differential equations use
       'rsolve' or 'dsolve' functions respectively, which are also
       wrappers combining set of problem specific methods.

       >>> from sympy import *
       >>> x, y, a = symbols('x', 'y', 'a')

       >>> solve(x**2 - 3*x + 2, x)
       [2, 1]

       >>> solve(x**2 == a, x)
       [-a**(1/2), a**(1/2)]

       #>>> solve(x**4 == 1, x) # use evalc() in polys
       #[-I, -1, 1, I]

       >>> solve([x + 5*y == 2, -3*x + 6*y == 15], [x, y])
       {y: 1, x: -3}

    """

    if isinstance(syms, Basic):
        syms = [syms]

    if not isinstance(eq, list):
        if isinstance(eq, Equality):
            # got equation, so move all the
            # terms to the left hand side
            equ = eq.lhs - eq.rhs
        else:
            equ = Basic.sympify(eq)

        try:
            # 'roots' method will return all possible complex
            # solutions, however we have to remove duplicates
            solutions = list(set(roots(equ, syms[0])))
        except PolynomialException:
            raise "Not a polynomial equation. Can't solve it, yet."

        if simplified == True:
            return [ simplify(s) for s in solutions ]
        else:
            return solutions
    else:
        # augmented matrix
        n, m = len(eq), len(syms)
        matrix = zeronm(n, m+1)

        index = {}

        for i in range(0, m):
            index[syms[i]] = i

        for i in range(0, n):
            if isinstance(eq[i], Equality):
                # got equation, so move all the
                # terms to the left hand side
                equ = eq[i].lhs - eq[i].rhs
            else:
                equ = Basic.sympify(eq[i])

            content = collect(equ, syms, evaluate=False)

            for var, expr in content.iteritems():
                if isinstance(var, Symbol) and not expr.has(*syms):
                    matrix[i, index[var]] = expr
                elif isinstance(var, Basic.One) and not expr.has(*syms):
                    matrix[i, m] = -expr
                else:
                    raise "Not a linear system. Can't solve it, yet."
        else:
            return solve_linear_system(matrix, syms, simplified)

def solve_linear_system(system, syms, simplified=True):
    """Solve system of N linear equations with M variables, which means
       both Cramer and over defined systems are supported. The possible
       number of solutions is zero, one or infinite. Respectively this
       functions will return empty dictionary or dictionary containing
       the same or less number of items as the number of variables. If
       there is infinite number of solutions, it will skip variables
       with can be assigned with arbitrary values.

       Input to this functions is a Nx(M+1) matrix, which means it has
       to be in augmented form. If you are unhappy with such setting
       use 'solve' method instead, where you can input equations
       explicitely. And don't worry aboute the matrix, this function
       is persistent and will make a local copy of it.

       The algorithm used here is fraction free Gaussian elimination,
       which results, after elimination, in upper-triangular matrix.
       Then solutions are found using back-substitution. This approach
       is more efficient and compact than the Gauss-Jordan method.

       >>> from sympy import *
       >>> x, y = symbols('x', 'y')

       Solve the following system:

              x + 4 y ==  2
           -2 x +   y == 14

       >>> system = Matrix(( (1, 4, 2), (-2, 1, 14)))
       >>> solve_linear_system(system, [x, y])
       {y: 2, x: -6}

    """
    matrix = system[:,:]     # we would like to be persistent
    i, m = 0, matrix.cols-1  # don't count augmentation

    while i < matrix.lines:
        if matrix [i, i] == 0:
            # there is no pivot in current column
            # so try to find one in other colums
            for k in range(i+1, m):
                if matrix[i, k] != 0:
                    break
            else:
                if matrix[i, m] != 0:
                    return {}   # no solutions
                else:
                    # zero row or was a linear combination of
                    # other rows so now we can safely skip it
                    matrix.row_del(i)
                    continue

            # we want to change the order of colums so
            # the order of variables must also change
            syms[i], syms[k] = syms[k], syms[i]
            matrix.col_swap(i, k)

        pivot = matrix [i, i]

        # divide all elements in the current row by the pivot
        matrix.row(i, lambda x, _: x / pivot)

        for k in range(i+1, matrix.lines):
            if matrix[k, i] != 0:
                coeff = matrix[k, i]

                # subtract from the current row the row containing
                # pivot and multiplied by extracted coefficient
                matrix.row(k, lambda x, j: x - matrix[i, j]*coeff)

        i += 1

    # if there weren't any problmes, augmented matrix is now
    # in row-echelon form so we can check how many solutions
    # there are and extract them using back substitution

    if len(syms) == matrix.lines:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in range(k+1, m):
                content -= matrix[k, j]*solutions[syms[j]]

            if simplified == True:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    elif len(syms) > matrix.lines:
        # this system will have infinite number of solutions
        # dependent on exactly len(syms) - i parameters
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in range(k+1, i):
                content -= matrix[k, j]*solutions[syms[j]]

            # run back-substitution for parameters
            for j in range(i, m):
                content -= matrix[k, j]*syms[j]

            if simplified == True:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    else:
        return {}   # no solutions

def solve_undetermined_coeffs(equ, coeffs, sym, simplified=True):
    """Solve equation of a type p(x; a_1, ..., a_k) == q(x) where both
       p, q are univariate polynomials and f depends on k parameters.
       The result of this functions is a dictionary with symbolic
       values of those parameters with respect to coefficiens in q.

       This functions accepts both Equations class instances and ordinary
       SymPy expressions. Specification of parameters and variable is
       obligatory for efficiency and simplicity reason.

       >>> from sympy import *
       >>> a, b, c, x = symbols('a', 'b', 'c', 'x')

       >>> solve_undetermined_coeffs(2*a*x + a+b == x, [a, b], x)
       {a: 1/2, b: -1/2}

       >>> solve_undetermined_coeffs(a*c*x + a+b == x, [a, b], x)
       {a: 1/c, b: -1/c}

    """
    if isinstance(equ, Equality):
        # got equation, so move all the
        # terms to the left hand side
        equ = equ.lhs - equ.rhs

    system = collect(equ, sym, evaluate=False).values()

    if not any([ equ.has(sym) for equ in system ]):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian ellimination algorithm
        return solve(system, coeffs, simplified)
    else:
        return {} # no solutions

def solve_linear_system_LU(matrix, syms):
    """ LU function works for invertible only """
    assert matrix.lines == matrix.cols-1
    A = matrix[:matrix.lines,:matrix.lines]
    b = matrix[:,matrix.cols-1:]
    soln = A.LUsolve(b)
    solutions = {}
    for i in range(soln.lines):
        solutions[syms[i]] = soln[i,0]
    return solutions
