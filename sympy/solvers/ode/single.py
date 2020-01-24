#
# This is the module for ODE solver classes for single ODEs.
#

from typing import ClassVar, Dict, List, Optional, Type

from sympy.core.expr import Expr
from sympy.core.function import AppliedUndef, Derivative, Function
from sympy.core.relational import Equality
from sympy.core.symbol import Symbol, Dummy
from sympy.integrals import Integral
from sympy.simplify import simplify

from sympy.solvers.solvers import solve


class ODEMatchError(NotImplementedError):
    """Raised if a SingleODESolver is asked to solve an ODE it does not match"""
    pass


class SingleODEProblem:
    """Represents an ordinary differential equation (ODE)

    This class is used internally in the by dsolve and related
    functions/classes so that properties of an ODE can be computed
    efficiently.

    Examples
    ========

    This class is used internally by dsolve. To instantiate an instance
    directly first define an ODE problem:

    >>> from sympy import Eq, Function, Symbol
    >>> x = Symbol('x')
    >>> f = Function('f')
    >>> eq = f(x).diff(x, 2)

    Now you can create a SingleODEProblem instance and query its properties:

    >>> from sympy.solvers.ode.single import SingleODEProblem
    >>> problem = SingleODEProblem(f(x).diff(x), f(x), x)
    >>> problem.eq
    Derivative(f(x), x)
    >>> problem.func
    f(x)
    >>> problem.sym
    x
    """

    # Instance attributes:
    eq = None  # type: Expr
    func = None  # type: AppliedUndef
    sym = None  # type: Symbol

    def __init__(self, eq, func, sym):
        self.eq = eq
        self.func = func
        self.sym = sym

    # TODO: Add methods that can be used by many ODE solvers:
    # order
    # is_linear()
    # get_linear_coefficients()
    # eq_prepared (the ODE in prepared form)


class SingleODESolver:
    """
    Base class for Single ODE solvers.

    Subclasses should implement the _matches and _get_general_solution
    methods. This class is not intended to be instantiated directly but its
    subclasses are as part of dsolve.

    Examples
    ========

    You can use a subclass of SingleODEProblem to solve a particular type of
    ODE. We first define a particular ODE problem:

    >>> from sympy import Eq, Function, Symbol
    >>> x = Symbol('x')
    >>> f = Function('f')
    >>> eq = f(x).diff(x, 2)

    Now we solve this problem using the NthAlgebraic solver which is a
    subclass of SingleODESolver:

    >>> from sympy.solvers.ode.single import NthAlgebraic, SingleODEProblem
    >>> problem = SingleODEProblem(eq, f(x), x)
    >>> solver = NthAlgebraic(problem)
    >>> solver.get_general_solution()
    [Eq(f(x), _C*x + _C)]

    The normal way to solve an ODE is to use dsolve (which would use
    NthAlgebraic and other solvers internally). When using dsolve a number of
    other things are done such as evaluating integrals, simplifying the
    solution and renumbering the constants:

    >>> from sympy import dsolve
    >>> dsolve(eq, hint='nth_algebraic')
    Eq(f(x), C1 + C2*x)
    """

    # Subclasses should store the hint name (the argument to dsolve) in this
    # attribute
    hint = None  # type: ClassVar[str]

    # Subclasses should define this to indicate if they support an _Integral
    # hint.
    has_integral = None  # type: ClassVar[bool]

    # The ODE to be solved
    ode_problem = None  # type: SingleODEProblem

    # Cache whether or not the equation has matched the method
    _matched = None  # type: Optional[bool]

    def __init__(self, ode_problem):
        self.ode_problem = ode_problem

    def matches(self) -> bool:
        if self._matched is None:
            self._matched = self._matches()
        return self._matched

    def get_general_solution(self, *, simplify: bool = True) -> List[Equality]:
        if not self.matches():
            msg = "%s solver can not solve:\n%s"
            raise ODEMatchError(msg % (self.hint, self.ode_problem.eq))
        return self._get_general_solution()

    def _matches(self) -> bool:
        msg = "Subclasses of SingleODESolver should implement matches."
        raise NotImplementedError(msg)

    def _get_general_solution(self, *, simplify: bool = True) -> List[Equality]:
        msg = "Subclasses of SingleODESolver should implement get_general_solution."
        raise NotImplementedError(msg)


class NthAlgebraic(SingleODESolver):
    r"""
    Solves an `n`\th order ordinary differential equation using algebra and
    integrals.

    There is no general form for the kind of equation that this can solve. The
    the equation is solved algebraically treating differentiation as an
    invertible algebraic function.

    Examples
    ========

    >>> from sympy import Function, dsolve, Eq
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> eq = Eq(f(x) * (f(x).diff(x)**2 - 1), 0)
    >>> dsolve(eq, f(x), hint='nth_algebraic')
    [Eq(f(x), 0), Eq(f(x), C1 - x), Eq(f(x), C1 + x)]

    Note that this solver can return algebraic solutions that do not have any
    integration constants (f(x) = 0 in the above example).
    """

    hint = 'nth_algebraic'
    has_integral = True  # nth_algebraic_Integral hint

    def _matches(self):
        r"""
        Matches any differential equation that nth_algebraic can solve. Uses
        `sympy.solve` but teaches it how to integrate derivatives.

        This involves calling `sympy.solve` and does most of the work of finding a
        solution (apart from evaluating the integrals).
        """
        eq = self.ode_problem.eq
        func = self.ode_problem.func
        var = self.ode_problem.sym

        # Derivative that solve can handle:
        diffx = self._get_diffx(var)

        # Replace derivatives wrt the independent variable with diffx
        def replace(eq, var):
            def expand_diffx(*args):
                differand, diffs = args[0], args[1:]
                toreplace = differand
                for v, n in diffs:
                    for _ in range(n):
                        if v == var:
                            toreplace = diffx(toreplace)
                        else:
                            toreplace = Derivative(toreplace, v)
                return toreplace
            return eq.replace(Derivative, expand_diffx)

        # Restore derivatives in solution afterwards
        def unreplace(eq, var):
            return eq.replace(diffx, lambda e: Derivative(e, var))

        subs_eqn = replace(eq, var)
        try:
            # turn off simplification to protect Integrals that have
            # _t instead of fx in them and would otherwise factor
            # as t_*Integral(1, x)
            solns = solve(subs_eqn, func, simplify=False)
        except NotImplementedError:
            solns = []

        solns = [simplify(unreplace(soln, var)) for soln in solns]
        solns = [Equality(func, soln) for soln in solns]

        self.solutions = solns
        return len(solns) != 0

    def _get_general_solution(self, *, simplify: bool = True):
        return self.solutions

    # This needs to produce an invertible function but the inverse depends
    # which variable we are integrating with respect to. Since the class can
    # be stored in cached results we need to ensure that we always get the
    # same class back for each particular integration variable so we store these
    # classes in a global dict:
    _diffx_stored = {}  # type: Dict[Symbol, Type[Function]]

    @staticmethod
    def _get_diffx(var):
        diffcls = NthAlgebraic._diffx_stored.get(var, None)

        if diffcls is None:
            # A class that behaves like Derivative wrt var but is "invertible".
            class diffx(Function):
                def inverse(self):
                    # don't use integrate here because fx has been replaced by _t
                    # in the equation; integrals will not be correct while solve
                    # is at work.
                    return lambda expr: Integral(expr, var) + Dummy('C')

            diffcls = NthAlgebraic._diffx_stored.setdefault(var, diffx)

        return diffcls
