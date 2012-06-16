import operator
from ode import dsolve
from solvers import solve
from sympy.core.containers import Tuple
from sympy.core.function import AppliedUndef, Derivative, Function, Lambda
from sympy.core.relational import Eq
from sympy.core.symbol import Symbol, IntConst
from sympy.matrices import Matrix
from sympy.utilities import numbered_symbols


def ode_system(exprs, funcs):
    """Solve a system of ODEs with optional initial conditions.

    For internal use only. The public interface is `dsolve`.

    It crashes or it just has undefined behavior in numerous cases. Check the
    TODOs in this docstring.

    Input arguments
    ===============

    - exprs: a list of expression equated to zero
    - funcs: a list of functions to solve for (should be already applied, i.e.
      `f(x)` and not `f`)

    Initial conditions
    ==================

    Initial conditions are supplied as part of the list of equations. For the
    moment the preprocessor for initial conditions is quite simple, thus take
    care to insure that there are no free variables in the equations
    corresponding to initial conditions.

    TODO
    Initial conditions for derivatives using Subs objects are not yet
    supported.

    Scope and Nonlinear equations
    =============================

    The coupled equations must be linear homogeneous. Equations that can be
    separated from the system can be of any type supported by the single
    equation solver in `dsolve`.

    Which equations can be separated is determined by checking which equations
    concern only a single function. This is simple and not very powerful
    method. It may be possible to extend it using the `solve` function.

    TODO
    For the moment this solver will probably crash on overdetermined
    nonconsistent systems.

    TODO
    It will also crash if the supplied coupled equation are nonlinear.

    TODO
    It will not work either for linear homogeneous systems that can not be
    diagonalized.

    TODO
    Nonhomogeneous linear system do not work either.

    Preprocessing
    =============

    Equations that can be separated (according to the abovementioned method
    used for separation) are not preprocessed and left to the single ODE
    solver. Equation that can not be separated are preprocessed by calling the
    `doit` method.

    Constants of integration
    ========================

    TODO
    As the single ode solve does not care that simplifying the constants can
    destroy information about the coupling we need to turn it off. This is not
    yet done.

    TODO
    Create smarter constant simplification routines.

    Implicit solutions
    ==================

    The single ODE solvers return Equality instances and not mappings for the
    solutions. This is necessary because sometimes an implicit solution
    (i.e. an algebraic equation in the function to be solved for) is easy to
    find and at the same time it is impossible to solve the resulting
    algebraic equation in order to obtain an explicit solution.

    TODO
    Currently this solver crashes if intermediate implicit solutions are
    generated.

    Output format
    =============

    For the sake of consistency with the single ODE solver, results are
    returned in the form of a list of Equality instances:

    - single function to solve for: list of Equations describing each
      existing solution.
    - multiple functions to solve for: list of systems of Equations describing
      each existing solution. The system of Equations is itself a list.

    TODO
    Currently the solver crashes if it needs to produce multiple solutions.

    Internal format
    ===============

    As implicit solutions are unsupported for the moment and as anyway
    supporting them for coupled equations is much more complicated than what we
    have here, there is no reason ot use Equality instances internally. Much
    simpler is to use dictionaries. They are transformend into lists of
    Equations at the end of the function.

    Examples
    ========
    Due to the different possible solutions these are not doctested.
    Check the test_ode_system.py file.
    """
    var = funcs[0].free_symbols.pop()
    init_conds   = [e for e in exprs if var not in e.free_symbols]
    separable_eq = [e for e in exprs if len(e.atoms(*funcs))==1
                                        and var in e.free_symbols]
    coupled_eq   = [e for e in exprs if len(e.atoms(*funcs))!=1
                                        and var in e.free_symbols]

    # Solve the separable equations. It is possible that after solving the
    # separable equations, some of the coupled equation become separable
    # themselves. This even takes care for triangular inhomogeneous systems.
    # However if the system is not already triangulized this leaves the work to
    # the next code block.
    separable_sols = {}
    while True:
        if separable_eq:
            sols = [dsolve(e) for e in separable_eq]
            separable_sols.update(dict((s.lhs, s.rhs) for s in sols))
        # Substitue the solved separable equations in the unsolved coupled ones.
        new_eq = [e.subs(separable_sols) for e in coupled_eq]
        # Check for changes.
        if new_eq == coupled_eq:
            break
        else:
            separable_eq = [e for e in new_eq if len(e.atoms(*funcs))==1]
            coupled_eq   = [e for e in new_eq if len(e.atoms(*funcs))!=1]

    # Solve the coupled equations.
    if coupled_eq:
        coupled_eq = [e.doit() for e in coupled_eq]
        coupled_funcs = set(funcs) - set(f for e in separable_eq
                                           for f in e.atoms(AppliedUndef))
        sols = ode_system_wo_ic(coupled_eq, list(coupled_funcs), var)
    else:
        sols = {}
    sols.update(separable_sols)

    # Solve the initial conditions.
    if init_conds:
        constants = list(set(c for e in sols.values()
                               for c in e.atoms(IntConst)))
        lambda_sols = dict((k.func, Lambda(var, v)) for k, v in sols.items())
        init_conds = [ic.subs(lambda_sols) for ic in init_conds]
        sub_init = solve(init_conds, constants)
        sols = dict((k, v.subs(sub_init)) for k, v in sols.items())

    # Not necessary for coding, however for aesthetic reasons renumerate the
    # constants.
    old_constants = list(set(c for e in sols.values()
                               for c in e.atoms(IntConst)))
    new_constants = [IntConst('C%d'%(i+1)) for i in range(len(old_constants))]
    sub_const = zip(old_constants, new_constants)
    sols = dict((k, v.subs(sub_const)) for k, v in sols.items())

    return [Eq(k, v) for k, v in sols.items()]


def ode_system_wo_ic(exprs, funcs, var):
    """Solve a linear homogeneous system of ODEs without initial conditions.

    The input `exprs` should not contain complicated derivative expressions.

    Returns a dictionary of the solutions. It is not necessary to return Eq
    instance because for linear homogeneous systems it is always possible to
    solve. Moreover this is for internal use, thus simplicity is more important
    than consistency with the public api.

    TODO: it fails on triangular systems
    TODO: it does not support nonhomogeneous systems

    >>> from sympy import Symbol, Function
    >>> from sympy.solvers.ode_system import ode_system_wo_ic
    >>> x = Symbol('x')
    >>> f = Function('f')(x)
    >>> g = Function('g')(x)
    >>> f_ = f.diff(x)
    >>> g_ = g.diff(x)
    >>> d = ode_system_wo_ic([f_+g, g_+f], [f, g], x)
    >>> d[g]-d[f]
    2*C1*exp(x)
    """
    exprs, funcs, subs = remove_higher_derivatives(exprs, funcs, var)
    matrix = construct_matrix(exprs, funcs, var)
    transf_matrix, diag_matrix = matrix.diagonalize()
    solutions = transf_matrix*ode_diagonal_system(diag_matrix, var)
    return dict([Tuple(*t).subs(subs) for t in zip(funcs, solutions)])


def ode_diagonal_system(diag_matrix, var):
    """Given a diagonal matrix, solve the coresponding ODEs.

    Returns the vector of solutions as expressions dependent on `var`.

    >>> from sympy import Symbol, Matrix
    >>> from sympy.solvers.ode_system import ode_diagonal_system
    >>> x = Symbol('x')
    >>> m = Matrix([[x, 0], [x**2, 0]])
    >>> ode_diagonal_system(m, x)
    [C1*exp(x**2/2)]
    [C1*exp(x**3/3)]
    """
    size = diag_matrix.shape[0]
    funcs = Matrix(size, 1, lambda a,b:dummy_func()(var))
    derivs = Matrix([f.diff(var) for f in funcs])
    equations = derivs - diag_matrix*funcs
    solutions = [dsolve(equ).rhs for equ in equations]
    return Matrix(solutions)


def construct_matrix(exprs, funcs, var):
    """For a system of first order ODEs, construct the coresponding matrix.

    The input `exprs` should not contain complicated derivative expressions.

    >>> from sympy import Symbol, Function
    >>> from sympy.solvers.ode_system import construct_matrix
    >>> x = Symbol('x')
    >>> f = Function('f')(x)
    >>> g = Function('g')(x)
    >>> f_ = f.diff(x)
    >>> g_ = g.diff(x)
    >>> construct_matrix([f_+g, g_+f], [f, g], x)
    [ 0, -1]
    [-1,  0]
    """
    derivs = [f.diff(var) for f in funcs]
    sol_dict = solve(exprs, derivs, dict=True)
    if sol_dict:
        sol_dict=sol_dict[0]
        matrix_content = [[sol_dict[d].coeff(f) for f in funcs]
                                                for d in derivs]
    else:
        raise ValueError('Inconsistent system of equations.')
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

    >>> from sympy import Symbol, Function
    >>> from sympy.solvers.ode_system import remove_higher_derivatives
    >>> x = Symbol('x')
    >>> f = Function('f')(x)
    >>> g = Function('g')(x)
    >>> f__ = f.diff(x, 2)
    >>> g_ = g.diff(x)
    >>> exprs, funcs, subs = remove_higher_derivatives([f__+g, g_+f], [f, g], x)
    >>> exprs
    [g(x) + Derivative(rrn1(x), x), f(x) + Derivative(g(x), x), rrn1(x) - Derivative(f(x), x)]
    >>> funcs
    [f(x), g(x), rrn1(x)]
    >>> subs
    {rrn1(x): Derivative(f(x), x)}
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
