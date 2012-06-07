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
    Initial conditions for derivatives using Subs objects are not yes
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

    >> from sympy import ode_system, Symbol, Function
    >> func = Function('f')
    >> gunc = Function('g')
    >> hunc = Function('h')
    >> x = Symbol('x')
    >> f = func(x)
    >> f_ = f.diff(x)
    >> f__ = f_.diff(x)
    >> g = gunc(x)
    >> g_ = g.diff(x)
    >> g__ = g_.diff(x)
    >> h = hunc(x)
    >> h_ = h.diff(x)
    >> h__ = h_.diff(x)

    Simple equation:
    >> sys = [f_+f]
    >> sol = ode_system(sys, [f])
    >> sol
    [f(x) == C1*exp(-x)]

    With initial conditions:
    >> sys = [f_+f, func(0)-2]
    >> sol = ode_system(sys, [f])
    >> sol
    [f(x) == 2*exp(-x)]

    With "boundary" conditions:
    >> sys = [f_+f, func(2)-3]
    >> sol = ode_system(sys, [f])
    >> sol
    [f(x) == 3*exp(2)*exp(-x)]

    Second order system:
    >> sys = [f__+f]
    >> sol = ode_system(sys, [f])
    >> sol
    [f(x) == C1*sin(x) + C2*cos(x)]

    Due to deficiencies in `Derivative` you must jump through hoops to get
    initials conditions working here (TODO use the Subs object):
    >> sys = [g-f_, f+g_] # The same as [f__+f] just substitute g_ = f__
    >> sol = ode_system(sys, [f, g])
    >> sol
    [g(x) == -I*C1*exp(-I*x) + I*C2*exp(I*x), f(x) == C1*exp(-I*x) + C2*exp(I*x)]

    Adding the initial conditions:
    >> sys = [g-f_, f+g_, func(0)-1, gunc(0)]
    >> sol = ode_system(sys, [f, g])
    >> sol
    [g(x) == I*exp(I*x)/2 - I*exp(-I*x)/2, f(x) == exp(I*x)/2 + exp(-I*x)/2]

    A separable nonlinear equation:
    >> sys = [f_+f**2, g_+h, h_-g]
    >> sol = ode_system(sys, [f,g,h])
    >> sol
    [g(x) == I*C1*exp(I*x) - I*C2*exp(-I*x), f(x) == 1/(C3 + x), h(x) == C1*exp(I*x) + C2*exp(-I*x)]

    A nonlinear system that is iteratively separable (we can separate and solve
    the first equation, substitute it in the second and solve it also, etc.):
    >> sys = [f_+f**2, g_*f-1]
    >> sol = ode_system(sys, [f,g])
    >> sol # Due to bad constant simplification there are too many constants
    [g(x) == C1 + C2*x + x**2/2, f(x) == 1/(C3 + x)]

    >> sys = [f_+f**2, g_*f-1, h_*f-g]
    >> sol = ode_system(sys, [f,g,h])
    >> sol # Due to bad constant simplification there are too many constants
    [g(x) == C4 + C5*x + x**2/2, f(x) == 1/(C1 + x), h(x) == C2*C3*x + C6 + x**4/8 + x**3*(C2/6 + C7/3) + x**2*(C2*C7/2 + C3/2)]

    """
    var = funcs[0].free_symbols.pop()
    init_conds   = [e for e in exprs if var not in e.free_symbols]
    separable_eq = [e for e in exprs if len(e.atoms(*funcs))==1
                                        and var in e.free_symbols]
    coupled_eq   = [e for e in exprs if len(e.atoms(*funcs))!=1
                                        and var in e.free_symbols]

    # Solve the separable equations. It is possible that after solving the
    # separable equations, some of the coupled equation become separable
    # themselves.
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

    # Not necessary for coding, however for aesthetic reason renumerate the
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
    than consistency with the public api."""
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
