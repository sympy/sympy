
""" This module contain solvers for all kinds of equations:

    - algebraic, use solve()

    - recurrence, use rsolve() (not implemented)

    - differential, use dsolve() (not implemented)

"""

from sympy.core.sympify import sympify
from sympy.core.basic import Basic, S, C
from sympy.core.symbol import Symbol, Wild
from sympy.core.relational import Equality
from sympy.core.function import Derivative, diff

from sympy.functions import sqrt, log, exp, LambertW
from sympy.simplify import simplify, collect
from sympy.matrices import Matrix, zeronm
from sympy.polys import roots

from sympy.utilities import any, all

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
       [I, 1, -1, -I]

       Solve a linear system:

       >>> solve((x+5*y-2, -3*x+6*y-15), x, y)
       {y: 1, x: -3}

    """
    if not symbols:
        raise ValueError('no symbols were given')

    if len(symbols) == 1:
        if isinstance(symbols[0], (list, tuple, set)):
            symbols = symbols[0]

    symbols = map(sympify, symbols)

    if any(not s.is_Symbol for s in symbols):
        raise TypeError('not a Symbol')

    if not isinstance(f, (tuple, list, set)):
        f = sympify(f)

        if isinstance(f, Equality):
            f = f.lhs - f.rhs

        if len(symbols) == 1:
            poly = f.as_poly(*symbols)

            if poly is not None:
                result = roots(poly, cubics=True, quartics=True).keys()
            else:
                result = [tsolve(f, *symbols)]
        else:
            raise NotImplementedError('multivariate equation')

        if flags.get('simplified', True):
            return map(simplify, result)
        else:
            return result
    else:
        if not f:
            return {}
        else:
            polys = []

            for g in f:
                g = sympify(g)

                if isinstance(g, Equality):
                    g = g.lhs - g.rhs

                poly = g.as_poly(*symbols)

                if poly is not None:
                    polys.append(poly)
                else:
                    raise NotImplementedError

            if all(p.is_linear for p in polys):
                n, m = len(f), len(symbols)
                matrix = zeronm(n, m + 1)

                for i, poly in enumerate(polys):
                    for coeff, monom in poly.iter_terms():
                        try:
                            j = list(monom).index(1)
                            matrix[i, j] = coeff
                        except ValueError:
                            matrix[i, m] = -coeff

                return solve_linear_system(matrix, *symbols, **flags)
            else:
                raise NotImplementedError('polynomial system')

def solve_linear_system(system, *symbols, **flags):
    """Solve system of N linear equations with M variables, which means
       both Cramer and over defined systems are supported. The possible
       number of solutions is zero, one or infinite. Respectively this
       procedure will return None or dictionary with solutions. In the
       case of over definend system all arbitrary parameters are skiped.
       This may cause situation in with empty dictionary is returned.
       In this case it means all symbols can be assigne arbitray values.

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
       >>> x, y = symbols('xy')

       Solve the following system:

              x + 4 y ==  2
           -2 x +   y == 14

       >>> system = Matrix(( (1, 4, 2), (-2, 1, 14)))
       >>> solve_linear_system(system, x, y)
       {y: 2, x: -6}

    """
    matrix = system[:,:]
    syms = list(symbols)

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

    simplified = flags.get('simplified', True)

    if len(syms) == matrix.lines:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in range(k+1, m):
                content -= matrix[k, j]*solutions[syms[j]]

            if simplified:
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
       values of those parameters with respect to coefficiens in q.

       This functions accepts both Equations class instances and ordinary
       SymPy expressions. Specification of parameters and variable is
       obligatory for efficiency and simplicity reason.

       >>> from sympy import *
       >>> a, b, c, x = symbols('a', 'b', 'c', 'x')

       >>> solve_undetermined_coeffs(Eq(2*a*x + a+b, x), [a, b], x)
       {b: -1/2, a: 1/2}

       >>> solve_undetermined_coeffs(Eq(a*c*x + a+b, x), [a, b], x)
       {b: -1/c, a: 1/c}

    """
    if isinstance(equ, Equality):
        # got equation, so move all the
        # terms to the left hand side
        equ = equ.lhs - equ.rhs

    system = collect(equ.expand(), sym, evaluate=False).values()

    if not any([ equ.has(sym) for equ in system ]):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian ellimination algorithm
        return solve(system, *coeffs, **flags)
    else:
        return None # no solutions

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

def dsolve(eq, funcs):
    """
    Solves any (supported) kind of differential equation.

    Usage
    =====
        dsolve(f, y(x)) -> Solve a differential equation f for the function y


    Details
    =======
        @param f: ordinary differential equation (either just the left hand
            side, or the Equality class)

        @param y: indeterminate function of one variable

        - you can declare the derivative of an unknown function this way:
        >>> from sympy import *
        >>> x = Symbol('x') # x is the independent variable

        >>> f = Function("f")(x) # f is a function of x
        >>> f_ = Derivative(f, x) # f_ will be the derivative of f with respect to x

        - This function just parses the equation "eq" and determines the type of
        differential equation by its order, then it determines all the coefficients and then
        calls the particular solver, which just accepts the coefficients.
        - "eq" can be either an Equality, or just the left hand side (in which
          case the right hand side is assumed to be 0)
        - see test_ode.py for many tests, that serve also as a set of examples
          how to use dsolve

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')

        >>> f = Function('f')
        >>> dsolve(Derivative(f(x),x,x)+9*f(x), f(x))
        C1*sin(3*x) + C2*cos(3*x)
        >>> dsolve(Eq(Derivative(f(x),x,x)+9*f(x)+1, 1), f(x))
        C1*sin(3*x) + C2*cos(3*x)

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

        x = f.args[0]
        f = f.func

        #We first get the order of the equation, so that we can choose the
        #corresponding methods. Currently, only first and second
        #order odes can be handled.
        order = deriv_degree(eq, f(x))

        if  order > 2 :
           raise NotImplementedError("dsolve: Cannot solve " + str(eq))
        elif order == 2:
            return solve_ODE_second_order(eq, f(x))
        elif order == 1:
            return solve_ODE_first_order(eq, f(x))
        else:
            raise NotImplementedError("Not a differential equation!")

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
    solves many kinds of first order odes, different methods are used
    depending on the form of the given equation. Now the linear
    case is implemented.
    """
    from sympy.integrals.integrals import integrate
    x = f.args[0]
    f = f.func

    #linear case: a(x)*f'(x)+b(x)*f(x)+c(x) = 0
    a = Wild('a', exclude=[f(x)])
    b = Wild('b', exclude=[f(x)])
    c = Wild('c', exclude=[f(x)])

    r = eq.match(a*diff(f(x),x) + b*f(x) + c)
    if r:
        t = C.exp(integrate(r[b]/r[a], x))
        tt = integrate(t*(-r[c]/r[a]), x)
        return (tt + Symbol("C1"))/t

    #other cases of first order odes will be implemented here

    raise NotImplementedError("dsolve: Cannot solve " + str(eq))

def solve_ODE_second_order(eq, f):
    """
    solves many kinds of second order odes, different methods are used
    depending on the form of the given equation. Now the constanst
    coefficients case and a special case are implemented.
    """
    x = f.args[0]
    f = f.func

    #constant coefficients case: af''(x)+bf'(x)+cf(x)=0
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])

    r = eq.match(a*f(x).diff(x,x) + c*f(x))
    if r:
        return Symbol("C1")*C.sin(sqrt(r[c]/r[a])*x)+Symbol("C2")*C.cos(sqrt(r[c]/r[a])*x)

    r = eq.match(a*f(x).diff(x,x) + b*diff(f(x),x) + c*f(x))
    if r:
        r1 = solve(r[a]*x**2 + r[b]*x + r[c], x)
        if r1[0].is_real:
            if len(r1) == 1:
                return (Symbol("C1") + Symbol("C2")*x)*exp(r1[0]*x)
            else:
                return Symbol("C1")*exp(r1[0]*x) + Symbol("C2")*exp(r1[1]*x)
        else:
            r2 = abs((r1[0] - r1[1])/(2*S.ImaginaryUnit))
            return (Symbol("C2")*C.cos(r2*x) + Symbol("C1")*C.sin(r2*x))*exp((r1[0] + r1[1])*x/2)

    #other cases of the second order odes will be implemented here

    #special equations, that we know how to solve
    t = x*C.exp(f(x))
    tt = a*t.diff(x, x)/t
    r = eq.match(tt.expand())
    if r:
        return -solve_ODE_1(f(x), x)

    t = x*C.exp(-f(x))
    tt = a*t.diff(x, x)/t
    r = eq.match(tt.expand())
    if r:
        #check, that we've rewritten the equation correctly:
        #assert ( r[a]*t.diff(x,2)/t ) == eq.subs(f, t)
        return solve_ODE_1(f(x), x)

    neq = eq*C.exp(f(x))/C.exp(-f(x))
    r = neq.match(tt.expand())
    if r:
        #check, that we've rewritten the equation correctly:
        #assert ( t.diff(x,2)*r[a]/t ).expand() == eq
        return solve_ODE_1(f(x), x)

    raise NotImplementedError("cannot solve this")

def solve_ODE_1(f, x):
    """ (x*exp(-f(x)))'' = 0 """
    C1 = Symbol("C1")
    C2 = Symbol("C2")
    return -C.log(C1+C2/x)

x = Symbol('x', dummy=True)
a,b,c,d,e,f,g,h = [Wild(t, exclude=[x]) for t in 'abcdefgh']
patterns = None

def _generate_patterns():
    """Generates patterns for transcendental equations.

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
        (1/2)/log(3)*(-5*log(3) + log(4))

        >>> tsolve(log(x) + 2*x, x)
        (1/2)*LambertW(2)

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
    # is finite), and we have a direct solution
    r = Wild('r')
    m = eq2.match((a*x+b)*r)
    if m and m[a]:
        return (-b/a).subs(m).subs(x, sym)
    for p, sol in patterns:
        m = eq2.match(p)
        if m:
            return sol.subs(m).subs(x, sym)
    raise ValueError("unable to solve the equation")
