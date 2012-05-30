import operator
from ode import dsolve
from solvers import solve
from sympy.core.containers import Tuple
from sympy.core.function import AppliedUndef, Derivative, Function, Lambda
from sympy.core.symbol import Symbol
from sympy.matrices import Matrix
from sympy.utilities import numbered_symbols


def ode_system(exprs, funcs, var):
    """Solve a system of ODEs with optional initial conditions.

    Initial conditions are supplied as part of the list of equations. For the
    moment the preprocessor for initial conditions is quite simple, thus take
    care to insure that there are no free variables in the equations
    corresponding to initial conditions.

    >>> from sympy import ode_system, Symbol, Function
    >>> func = Function('f')
    >>> gunc = Function('g')
    >>> x = Symbol('x')
    >>> f = func(x)
    >>> f_ = f.diff(x)
    >>> f__ = f_.diff(x)
    >>> g = gunc(x)
    >>> g_ = g.diff(x)
    >>> g__ = g_.diff(x)

    Simple equation:
    >>> sys = [f_+f]
    >>> sol = ode_system(sys, [f], x)
    >>> sol[f]
    C1*exp(-x)

    With initial conditions:
    >>> sys = [f_+f, func(0)-2]
    >>> sol = ode_system(sys, [f], x)
    >>> sol[f]
    2*exp(-x)

    With "boundary" conditions:
    >>> sys = [f_+f, func(2)-3]
    >>> sol = ode_system(sys, [f], x)
    >>> sol[f]
    3*exp(2)*exp(-x)

    Second order system:
    >>> sys = [f__+f]
    >>> sol = ode_system(sys, [f], x)
    >>> sol[f]
    I*C1*exp(-I*x) - I*C2*exp(I*x)

    Due to deficiencies in `Derivative` you must jump through hoops to get
    initials conditions working here:
    >>> sys = [g-f_, f+g_] # The same as [f__+f] just substitute g_ = f__
    >>> sol = ode_system(sys, [f, g], x)
    >>> sol[f]
    I*C1*exp(-I*x) - I*C2*exp(I*x)

    Adding the initial conditions:
    >>> sys = [g-f_, f+g_, func(0)-1, gunc(0)]
    >>> sol = ode_system(sys, [f, g], x)
    >>> sol[f]
    exp(I*x)/2 + exp(-I*x)/2
    """
    # TODO preprocessor the `doit` the derivatives

    init_conds = [e for e in exprs if var not in e.free_symbols]
    diff_eq = [e for e in exprs if var in e.free_symbols]
    sols = ode_system_wo_ic(diff_eq, funcs, var)
    # XXX TODO Workaround: it is difficult to collect the constants
    constants = list(set(c for sol in sols.values()
                           for c in sol.free_symbols
                           if c.name.startswith('C')))
    # XXX End of workaround.
    lambda_sols = dict((k.func, Lambda(var, v)) for k, v in sols.items())
    init_conds = [ic.subs(lambda_sols) for ic in init_conds]
    c = solve(init_conds, constants)
    sols = dict((k, v.subs(c)) for k, v in sols.items())
    return sols


def ode_system_wo_ic(exprs, funcs, var):
    """Solve a system of ODEs without initial conditions.

    The input `exprs` should not contain complicated derivative expressions.

    Returns a dictionary of the solutions."""
    exprs, funcs, subs = remove_higher_derivatives(exprs, funcs, var)
    matrix = construct_matrix(exprs, funcs, var)
    transf_matrix, diag_matrix = matrix.diagonalize()
    solutions = transf_matrix*ode_diagonal_system(diag_matrix, var)
    return dict([Tuple(*t).subs(subs) for t in zip(funcs, solutions)])


def ode_diagonal_system(diag_matrix, var):
    """Given a diagonal matrix, solve the coresponding ODEs.

    Returns the vector of solutions as expressions dependent on `var`."""
    size = diag_matrix.shape[0]
    funcs = Matrix(size, 1, lambda a,b:dummy_func()(var))
    derivs = Matrix([f.diff(var) for f in funcs])
    equations = derivs - diag_matrix*funcs
    solutions = [dsolve(equ).rhs for equ in equations]
    # XXX TODO Workaround: new independent constants
    counter = 1
    for i, sol in enumerate(solutions):
        current_constants = [c for c in sol.free_symbols if c.name.startswith('C')]
        solutions[i] = sol.subs(zip(current_constants, numbered_symbols(prefix='C', cls=Symbol, start=counter)))
        counter += len(current_constants)
    # XXX End of workaround.
    return Matrix(solutions)


def construct_matrix(exprs, funcs, var):
    """For a system of first order ODEs, construct the coresponding matrix.

    The input `exprs` should not contain complicated derivative expressions.
    """
    derivs = [f.diff(var) for f in funcs]
    sol_dict = solve(exprs, derivs)
    matrix_content = [[sol_dict[d].coeff(f) for f in funcs]
                                            for d in derivs]
    return Matrix(matrix_content)


def remove_higher_derivatives(exprs, funcs, var, subs={}):
    """Substitute higher derivatives with first derivatives of dummy functions.

    Used as preprocessor for ode_system. The input `exprs` should not contain
    complicated derivative expressions.

    Argumets:
    =========

    - exprs - the equations in the form of expressions assumed equal to zero
    - funcs - the function wrt which to solve
    - var   - the independent variable
    - subs  - substitutions to be added to the return value (needed for recursion)

    Returns:
    ========

    - exprs - the higher derivatives were substituted and additional equation were added
    - funcs - the functions wrt which to solve with the new dummy function appended at the end
    - subs  - dictionarry of substitutions to remove the dummy functions
    """
    higher_order = list(set(d for e in exprs for d in e.atoms(Derivative)
                              if order_wrt(d, var)>1))
    if not higher_order:
        return exprs, funcs, subs
    else:
        one_lesser_order = [d.integrate(var) for d in higher_order]

        new_funcs = [dummy_func()(var) for d in higher_order]
        new_funcs_diff = [f.diff(var) for f in new_funcs]

        new_expr = [operator.sub(*t) for t in zip(new_funcs, one_lesser_order)]
        higher_to_lower = zip(higher_order, new_funcs_diff)
        mod_expr = [e.subs(higher_to_lower) for e in exprs]

        i_subs = dict(zip(new_funcs, one_lesser_order))
        i_subs.update(subs)

        return remove_higher_derivatives(mod_expr + new_expr,
                                         funcs + new_funcs,
                                         var,
                                         i_subs)


def order_wrt(derivative, symbol):
    """Return the order of a derivative wrt a symbol.

    The derivative should not be complicated."""
    return sum(symbol==a for a in derivative.args)


def dummy_func(counter=0):
    """Yeah, sure... This will blow up in my face."""
    # TODO Needs a real Dummy function.
    # XXX rrn stands for 'Really Random Name'. It is an ISO standard.
    counter += 1
    return Function('rrn%d' % counter)
