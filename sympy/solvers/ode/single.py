#
# This is the module for ODE solver classes for single ODEs.
#

import typing
from collections import defaultdict
if typing.TYPE_CHECKING:
    from typing import ClassVar
from typing import Dict, Type

from typing import Iterator, List, Optional
from sympy.core import Add, S, Pow
from sympy.core.exprtools import factor_terms
from sympy.core.expr import Expr
from sympy.core.function import AppliedUndef, Derivative, Function, expand, Subs, _mexpand, expand_mul
from sympy.core.numbers import Float, zoo
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Symbol, Dummy, Wild
from sympy.core.mul import Mul
from sympy.functions import exp, tan, cos, cosh, im, log, re, sin, sinh, sqrt, \
    atan2, conjugate
from sympy.integrals import Integral
from sympy.polys import (Poly, rootof, roots)
from sympy.polys.polytools import cancel, factor
from sympy.simplify import collect, simplify, separatevars, logcombine, powsimp, trigsimp
from sympy.simplify.radsimp import fraction
from sympy.utilities import numbered_symbols, default_sort_key
from sympy.solvers.solvers import solve
from sympy.matrices import wronskian
from sympy.solvers.deutils import ode_order, _preprocess
from .subscheck import sub_func_doit
from .hypergeometric import equivalence_hypergeometric, match_2nd_2F1_hypergeometric, get_sol_2F1_hypergeometric, match_2nd_hypergeometric


class ODEMatchError(NotImplementedError):
    """Raised if a SingleODESolver is asked to solve an ODE it does not match"""
    pass


def cached_property(func):
    '''Decorator to cache property method'''
    attrname = '_' + func.__name__
    def propfunc(self):
        val = getattr(self, attrname, None)
        if val is None:
            val = func(self)
            setattr(self, attrname, val)
        return val
    return property(propfunc)


class SingleODEProblem:
    """Represents an ordinary differential equation (ODE)

    This class is used internally in the by dsolve and related
    functions/classes so that properties of an ODE can be computed
    efficiently.

    Examples
    ========

    This class is used internally by dsolve. To instantiate an instance
    directly first define an ODE problem:

    >>> from sympy import Function, Symbol
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
    _order = None  # type: int
    _eq_expanded = None  # type: Expr
    _eq_preprocessed = None  # type: Expr
    _eq_high_order_free = None

    def __init__(self, eq, func, sym, prep=True):
        assert isinstance(eq, Expr)
        assert isinstance(func, AppliedUndef)
        assert isinstance(sym, Symbol)
        assert isinstance(prep, bool)
        self.eq = eq
        self.func = func
        self.sym = sym
        self.prep = prep

    @cached_property
    def order(self) -> int:
        return ode_order(self.eq, self.func)

    @cached_property
    def eq_preprocessed(self) -> Expr:
        return self._get_eq_preprocessed()

    @cached_property
    def eq_high_order_free(self) -> Expr:
        a = Wild('a', exclude=[self.func])
        c1 = Wild('c1', exclude=[self.sym])
        self.eq = expand(self.eq)
        # Precondition to try remove f(x) from highest order derivative
        reduced_eq = None
        if self.eq.is_Add:
            deriv_coef = self.eq.coeff(self.func.diff(self.sym, self.order))
            if deriv_coef not in (1, 0):
                r = deriv_coef.match(a*self.func**c1)
                if r and r[c1]:
                    den = self.func**r[c1]
                    reduced_eq = Add(*[arg/den for arg in self.eq.args])
        if not reduced_eq:
            reduced_eq = self.eq
        return reduced_eq

    @cached_property
    def eq_expanded(self) -> Expr:
        return expand(self.eq_preprocessed)

    def _get_eq_preprocessed(self) -> Expr:
        if self.prep:
            process_eq, process_func = _preprocess(self.eq, self.func)
            if process_func != self.func:
                raise ValueError
        else:
            process_eq = self.eq
        return process_eq

    def get_numbered_constants(self, num=1, start=1, prefix='C') -> List[Symbol]:
        """
        Returns a list of constants that do not occur
        in eq already.
        """
        ncs = self.iter_numbered_constants(start, prefix)
        Cs = [next(ncs) for i in range(num)]
        return Cs

    def iter_numbered_constants(self, start=1, prefix='C') -> Iterator[Symbol]:
        """
        Returns an iterator of constants that do not occur
        in eq already.
        """
        atom_set = self.eq.free_symbols
        func_set = self.eq.atoms(Function)
        if func_set:
            atom_set |= {Symbol(str(f.func)) for f in func_set}
        return numbered_symbols(start=start, prefix=prefix, exclude=atom_set)

    @cached_property
    def is_autonomous(self):
        u = Dummy('u')
        x = self.sym
        syms = self.eq.subs(self.func, u).free_symbols
        return x not in syms

    def get_linear_coefficients(self, eq, func, order):
        r"""
        Matches a differential equation to the linear form:

        .. math:: a_n(x) y^{(n)} + \cdots + a_1(x)y' + a_0(x) y + B(x) = 0

        Returns a dict of order:coeff terms, where order is the order of the
        derivative on each term, and coeff is the coefficient of that derivative.
        The key ``-1`` holds the function `B(x)`. Returns ``None`` if the ODE is
        not linear.  This function assumes that ``func`` has already been checked
        to be good.

        Examples
        ========

        >>> from sympy import Function, cos, sin
        >>> from sympy.abc import x
        >>> from sympy.solvers.ode.ode import _nth_linear_match
        >>> f = Function('f')
        >>> _nth_linear_match(f(x).diff(x, 3) + 2*f(x).diff(x) +
        ... x*f(x).diff(x, 2) + cos(x)*f(x).diff(x) + x - f(x) -
        ... sin(x), f(x), 3)
        {-1: x - sin(x), 0: -1, 1: cos(x) + 2, 2: x, 3: 1}
        >>> _nth_linear_match(f(x).diff(x, 3) + 2*f(x).diff(x) +
        ... x*f(x).diff(x, 2) + cos(x)*f(x).diff(x) + x - f(x) -
        ... sin(f(x)), f(x), 3) == None
        True

        """
        x = func.args[0]
        one_x = {x}
        terms = {i: S.Zero for i in range(-1, order + 1)}
        for i in Add.make_args(eq):
            if not i.has(func):
                terms[-1] += i
            else:
                c, f = i.as_independent(func)
                if (isinstance(f, Derivative)
                        and set(f.variables) == one_x
                        and f.args[0] == func):
                    terms[f.derivative_count] += c
                elif f == func:
                    terms[len(f.args[1:])] += c
                else:
                    return None
        return terms

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

    >>> from sympy import Function, Symbol
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

    # Subclasses should store in this attribute the list of order(s) of ODE
    # that subclass can solve or leave it to None if not specific to any order
    order = None  # type: Optional[list]

    def __init__(self, ode_problem):
        self.ode_problem = ode_problem

    def matches(self) -> bool:
        if self.order is not None and self.ode_problem.order not in self.order:
            self._matched = False
            return self._matched

        if self._matched is None:
            self._matched = self._matches()
        return self._matched

    def get_general_solution(self, *, simplify: bool = True) -> List[Equality]:
        if not self.matches():
            msg = "%s solver can not solve:\n%s"
            raise ODEMatchError(msg % (self.hint, self.ode_problem.eq))
        return self._get_general_solution(simplify_flag=simplify)

    def _matches(self) -> bool:
        msg = "Subclasses of SingleODESolver should implement matches."
        raise NotImplementedError(msg)

    def _get_general_solution(self, *, simplify_flag: bool = True) -> List[Equality]:
        msg = "Subclasses of SingleODESolver should implement get_general_solution."
        raise NotImplementedError(msg)


class SinglePatternODESolver(SingleODESolver):
    '''Superclass for ODE solvers based on pattern matching'''

    def wilds(self):
        prob = self.ode_problem
        f = prob.func.func
        x = prob.sym
        order = prob.order
        return self._wilds(f, x, order)

    def wilds_match(self):
        match = self._wilds_match
        return [match.get(w, S.Zero) for w in self.wilds()]

    def _matches(self):
        eq = self.ode_problem.eq_expanded
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        order = self.ode_problem.order
        df = f(x).diff(x, order)

        if order not in [1, 2]:
            return False

        pattern = self._equation(f(x), x, order)

        if not pattern.coeff(df).has(Wild):
            eq = expand(eq / eq.coeff(df))
        eq = eq.collect([f(x).diff(x), f(x)], func = cancel)

        self._wilds_match = match = eq.match(pattern)
        if match is not None:
            return self._verify(f(x))
        return False

    def _verify(self, fx) -> bool:
        return True

    def _wilds(self, f, x, order):
        msg = "Subclasses of SingleODESolver should implement _wilds"
        raise NotImplementedError(msg)

    def _equation(self, fx, x, order):
        msg = "Subclasses of SingleODESolver should implement _equation"
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

    def _get_general_solution(self, *, simplify_flag: bool = True):
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


class FirstExact(SinglePatternODESolver):
    r"""
    Solves 1st order exact ordinary differential equations.

    A 1st order differential equation is called exact if it is the total
    differential of a function. That is, the differential equation

    .. math:: P(x, y) \,\partial{}x + Q(x, y) \,\partial{}y = 0

    is exact if there is some function `F(x, y)` such that `P(x, y) =
    \partial{}F/\partial{}x` and `Q(x, y) = \partial{}F/\partial{}y`.  It can
    be shown that a necessary and sufficient condition for a first order ODE
    to be exact is that `\partial{}P/\partial{}y = \partial{}Q/\partial{}x`.
    Then, the solution will be as given below::

        >>> from sympy import Function, Eq, Integral, symbols, pprint
        >>> x, y, t, x0, y0, C1= symbols('x,y,t,x0,y0,C1')
        >>> P, Q, F= map(Function, ['P', 'Q', 'F'])
        >>> pprint(Eq(Eq(F(x, y), Integral(P(t, y), (t, x0, x)) +
        ... Integral(Q(x0, t), (t, y0, y))), C1))
                    x                y
                    /                /
                   |                |
        F(x, y) =  |  P(t, y) dt +  |  Q(x0, t) dt = C1
                   |                |
                  /                /
                  x0               y0

    Where the first partials of `P` and `Q` exist and are continuous in a
    simply connected region.

    A note: SymPy currently has no way to represent inert substitution on an
    expression, so the hint ``1st_exact_Integral`` will return an integral
    with `dy`.  This is supposed to represent the function that you are
    solving for.

    Examples
    ========

    >>> from sympy import Function, dsolve, cos, sin
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x),
    ... f(x), hint='1st_exact')
    Eq(x*cos(f(x)) + f(x)**3/3, C1)

    References
    ==========

    - https://en.wikipedia.org/wiki/Exact_differential_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 73

    # indirect doctest

    """
    hint = "1st_exact"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        P = Wild('P', exclude=[f(x).diff(x)])
        Q = Wild('Q', exclude=[f(x).diff(x)])
        return P, Q

    def _equation(self, fx, x, order):
        P, Q = self.wilds()
        return P + Q*fx.diff(x)

    def _verify(self, fx) -> bool:
        P, Q = self.wilds()
        x = self.ode_problem.sym
        y = Dummy('y')

        m, n = self.wilds_match()

        m = m.subs(fx, y)
        n = n.subs(fx, y)
        numerator = cancel(m.diff(y) - n.diff(x))

        if numerator.is_zero:
            # Is exact
            return True
        else:
            # The following few conditions try to convert a non-exact
            # differential equation into an exact one.
            # References:
            # 1. Differential equations with applications
            # and historical notes - George E. Simmons
            # 2. https://math.okstate.edu/people/binegar/2233-S99/2233-l12.pdf

            factor_n = cancel(numerator/n)
            factor_m = cancel(-numerator/m)
            if y not in factor_n.free_symbols:
                # If (dP/dy - dQ/dx) / Q = f(x)
                # then exp(integral(f(x))*equation becomes exact
                factor = factor_n
                integration_variable = x
            elif x not in factor_m.free_symbols:
                # If (dP/dy - dQ/dx) / -P = f(y)
                # then exp(integral(f(y))*equation becomes exact
                factor = factor_m
                integration_variable = y
            else:
                # Couldn't convert to exact
                return False

            factor = exp(Integral(factor, integration_variable))
            m *= factor
            n *= factor
            self._wilds_match[P] = m.subs(y, fx)
            self._wilds_match[Q] = n.subs(y, fx)
            return True

    def _get_general_solution(self, *, simplify_flag: bool = True):
        m, n = self.wilds_match()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        (C1,) = self.ode_problem.get_numbered_constants(num=1)
        y = Dummy('y')

        m = m.subs(fx, y)
        n = n.subs(fx, y)

        gen_sol = Eq(Subs(Integral(m, x)
                          + Integral(n - Integral(m, x).diff(y), y), y, fx), C1)
        return [gen_sol]


class FirstLinear(SinglePatternODESolver):
    r"""
    Solves 1st order linear differential equations.

    These are differential equations of the form

    .. math:: dy/dx + P(x) y = Q(x)\text{.}

    These kinds of differential equations can be solved in a general way.  The
    integrating factor `e^{\int P(x) \,dx}` will turn the equation into a
    separable equation.  The general solution is::

        >>> from sympy import Function, dsolve, Eq, pprint, diff, sin
        >>> from sympy.abc import x
        >>> f, P, Q = map(Function, ['f', 'P', 'Q'])
        >>> genform = Eq(f(x).diff(x) + P(x)*f(x), Q(x))
        >>> pprint(genform)
                    d
        P(x)*f(x) + --(f(x)) = Q(x)
                    dx
        >>> pprint(dsolve(genform, f(x), hint='1st_linear_Integral'))
                /       /                   \
                |      |                    |
                |      |         /          |     /
                |      |        |           |    |
                |      |        | P(x) dx   |  - | P(x) dx
                |      |        |           |    |
                |      |       /            |   /
        f(x) = |C1 +  | Q(x)*e           dx|*e
                |      |                    |
                \     /                     /


    Examples
    ========

    >>> f = Function('f')
    >>> pprint(dsolve(Eq(x*diff(f(x), x) - f(x), x**2*sin(x)),
    ... f(x), '1st_linear'))
    f(x) = x*(C1 - cos(x))

    References
    ==========

    - https://en.wikipedia.org/wiki/Linear_differential_equation#First_order_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 92

    # indirect doctest

    """
    hint = '1st_linear'
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        P = Wild('P', exclude=[f(x)])
        Q = Wild('Q', exclude=[f(x), f(x).diff(x)])
        return P, Q

    def _equation(self, fx, x, order):
        P, Q = self.wilds()
        return fx.diff(x) + P*fx - Q

    def _get_general_solution(self, *, simplify_flag: bool = True):
        P, Q = self.wilds_match()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        (C1,)  = self.ode_problem.get_numbered_constants(num=1)
        gensol = Eq(fx, ((C1 + Integral(Q*exp(Integral(P, x)),x))
            * exp(-Integral(P, x))))
        return [gensol]


class AlmostLinear(SinglePatternODESolver):
    r"""
    Solves an almost-linear differential equation.

    The general form of an almost linear differential equation is

    .. math:: a(x) g'(f(x)) f'(x) + b(x) g(f(x)) + c(x)

    Here `f(x)` is the function to be solved for (the dependent variable).
    The substitution `g(f(x)) = u(x)` leads to a linear differential equation
    for `u(x)` of the form `a(x) u' + b(x) u + c(x) = 0`. This can be solved
    for `u(x)` by the `first_linear` hint and then `f(x)` is found by solving
    `g(f(x)) = u(x)`.

    See Also
    ========
    :obj:`sympy.solvers.ode.single.FirstLinear`

    Examples
    ========

    >>> from sympy import Function, pprint, sin, cos
    >>> from sympy.solvers.ode import dsolve
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> d = f(x).diff(x)
    >>> eq = x*d + x*f(x) + 1
    >>> dsolve(eq, f(x), hint='almost_linear')
    Eq(f(x), (C1 - Ei(x))*exp(-x))
    >>> pprint(dsolve(eq, f(x), hint='almost_linear'))
                        -x
    f(x) = (C1 - Ei(x))*e
    >>> example = cos(f(x))*f(x).diff(x) + sin(f(x)) + 1
    >>> pprint(example)
                        d
    sin(f(x)) + cos(f(x))*--(f(x)) + 1
                        dx
    >>> pprint(dsolve(example, f(x), hint='almost_linear'))
                    /    -x    \             /    -x    \
    [f(x) = pi - asin\C1*e   - 1/, f(x) = asin\C1*e   - 1/]


    References
    ==========

    - Joel Moses, "Symbolic Integration - The Stormy Decade", Communications
      of the ACM, Volume 14, Number 8, August 1971, pp. 558
    """
    hint = "almost_linear"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        P = Wild('P', exclude=[f(x).diff(x)])
        Q = Wild('Q', exclude=[f(x).diff(x)])
        return P, Q

    def _equation(self, fx, x, order):
        P, Q = self.wilds()
        return P*fx.diff(x) + Q

    def _verify(self, fx):
        a, b = self.wilds_match()
        c, b = b.as_independent(fx) if b.is_Add else (S.Zero, b)
        # a, b and c are the function a(x), b(x) and c(x) respectively.
        # c(x) is obtained by separating out b as terms with and without fx i.e, l(y)
        # The following conditions checks if the given equation is an almost-linear differential equation using the fact that
        # a(x)*(l(y))' / l(y)' is independent of l(y)

        if b.diff(fx) != 0 and not simplify(b.diff(fx)/a).has(fx):
            self.ly = factor_terms(b).as_independent(fx, as_Add=False)[1] # Gives the term containing fx i.e., l(y)
            self.ax = a / self.ly.diff(fx)
            self.cx = -c  # cx is taken as -c(x) to simplify expression in the solution integral
            self.bx = factor_terms(b) / self.ly
            return True

        return False

    def _get_general_solution(self, *, simplify_flag: bool = True):
        x = self.ode_problem.sym
        (C1,)  = self.ode_problem.get_numbered_constants(num=1)
        gensol = Eq(self.ly, ((C1 + Integral((self.cx/self.ax)*exp(Integral(self.bx/self.ax, x)),x))
                * exp(-Integral(self.bx/self.ax, x))))

        return [gensol]


class Bernoulli(SinglePatternODESolver):
    r"""
    Solves Bernoulli differential equations.

    These are equations of the form

    .. math:: dy/dx + P(x) y = Q(x) y^n\text{, }n \ne 1`\text{.}

    The substitution `w = 1/y^{1-n}` will transform an equation of this form
    into one that is linear (see the docstring of
    :obj:`~sympy.solvers.ode.single.FirstLinear`).  The general solution is::

        >>> from sympy import Function, dsolve, Eq, pprint
        >>> from sympy.abc import x, n
        >>> f, P, Q = map(Function, ['f', 'P', 'Q'])
        >>> genform = Eq(f(x).diff(x) + P(x)*f(x), Q(x)*f(x)**n)
        >>> pprint(genform)
                    d                n
        P(x)*f(x) + --(f(x)) = Q(x)*f (x)
                    dx
        >>> pprint(dsolve(genform, f(x), hint='Bernoulli_Integral'), num_columns=110)
                                                                                                              -1
                                                                                                             -----
                                                                                                             n - 1
               //         /                                /                           \                    \
               ||        |                                |                            |                    |
               ||        |                 /              |                 /          |            /       |
               ||        |                |               |                |           |           |        |
               ||        |       (1 - n)* | P(x) dx       |       (1 - n)* | P(x) dx   |  (n - 1)* | P(x) dx|
               ||        |                |               |                |           |           |        |
               ||        |               /                |               /            |          /         |
        f(x) = ||C1 - n* | Q(x)*e                   dx +  | Q(x)*e                   dx|*e                  |
               ||        |                                |                            |                    |
               \\       /                                /                             /                    /


    Note that the equation is separable when `n = 1` (see the docstring of
    :obj:`~sympy.solvers.ode.single.Separable`).

    >>> pprint(dsolve(Eq(f(x).diff(x) + P(x)*f(x), Q(x)*f(x)), f(x),
    ... hint='separable_Integral'))
    f(x)
        /
    |                /
    |  1            |
    |  - dy = C1 +  | (-P(x) + Q(x)) dx
    |  y            |
    |              /
    /


    Examples
    ========

    >>> from sympy import Function, dsolve, Eq, pprint, log
    >>> from sympy.abc import x
    >>> f = Function('f')

    >>> pprint(dsolve(Eq(x*f(x).diff(x) + f(x), log(x)*f(x)**2),
    ... f(x), hint='Bernoulli'))
                    1
    f(x) =  -----------------
            C1*x + log(x) + 1

    References
    ==========

    - https://en.wikipedia.org/wiki/Bernoulli_differential_equation

    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 95

    # indirect doctest

    """
    hint = "Bernoulli"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        P = Wild('P', exclude=[f(x)])
        Q = Wild('Q', exclude=[f(x)])
        n = Wild('n', exclude=[x, f(x), f(x).diff(x)])
        return P, Q, n

    def _equation(self, fx, x, order):
        P, Q, n = self.wilds()
        return fx.diff(x) + P*fx - Q*fx**n

    def _get_general_solution(self, *, simplify_flag: bool = True):
        P, Q, n = self.wilds_match()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        (C1,) = self.ode_problem.get_numbered_constants(num=1)
        if n==1:
            gensol = Eq(log(fx), (
            C1 + Integral((-P + Q),x)
        ))
        else:
            gensol = Eq(fx**(1-n), (
                (C1 - (n - 1) * Integral(Q*exp(-n*Integral(P, x))
                            * exp(Integral(P, x)), x)
                ) * exp(-(1 - n)*Integral(P, x)))
            )
        return [gensol]


class Factorable(SingleODESolver):
    r"""
        Solves equations having a solvable factor.

        This function is used to solve the equation having factors. Factors may be of type algebraic or ode. It
        will try to solve each factor independently. Factors will be solved by calling dsolve. We will return the
        list of solutions.

        Examples
        ========

        >>> from sympy import Function, dsolve, pprint
        >>> from sympy.abc import x
        >>> f = Function('f')
        >>> eq = (f(x)**2-4)*(f(x).diff(x)+f(x))
        >>> pprint(dsolve(eq, f(x)))
                                        -x
        [f(x) = 2, f(x) = -2, f(x) = C1*e  ]


        """
    hint = "factorable"
    has_integral = False

    def _matches(self):
        eq = self.ode_problem.eq
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        order =self.ode_problem.order
        df = f(x).diff(x)
        self.eqs = []
        eq = eq.collect(f(x), func = cancel)
        eq = fraction(factor(eq))[0]
        factors = Mul.make_args(factor(eq))
        roots = [fac.as_base_exp() for fac in factors if len(fac.args)!=0]
        if len(roots)>1 or roots[0][1]>1:
            for base,expo in roots:
                if base.has(f(x)):
                    self.eqs.append(base)
            if len(self.eqs)>0:
                return True
        roots = solve(eq, df)
        if len(roots)>0:
            self.eqs = [(df - root) for root in roots]
            if len(self.eqs)==1:
                if order>1:
                    return False
                if self.eqs[0].has(Float):
                    return False
                return fraction(factor(self.eqs[0]))[0]-eq!=0
            return True
        return False


    def _get_general_solution(self, *, simplify_flag: bool = True):
        func = self.ode_problem.func.func
        x = self.ode_problem.sym
        eqns = self.eqs
        sols = []
        for eq in eqns:
            try:
                sol = dsolve(eq, func(x))
            except NotImplementedError:
                continue
            else:
                if isinstance(sol, list):
                    sols.extend(sol)
                else:
                    sols.append(sol)

        if sols == []:
            raise NotImplementedError("The given ODE " + str(eq) + " cannot be solved by"
                + " the factorable group method")
        return sols


class RiccatiSpecial(SinglePatternODESolver):
    r"""
    The general Riccati equation has the form

    .. math:: dy/dx = f(x) y^2 + g(x) y + h(x)\text{.}

    While it does not have a general solution [1], the "special" form, `dy/dx
    = a y^2 - b x^c`, does have solutions in many cases [2].  This routine
    returns a solution for `a(dy/dx) = b y^2 + c y/x + d/x^2` that is obtained
    by using a suitable change of variables to reduce it to the special form
    and is valid when neither `a` nor `b` are zero and either `c` or `d` is
    zero.

    >>> from sympy.abc import x, a, b, c, d
    >>> from sympy.solvers.ode import dsolve, checkodesol
    >>> from sympy import pprint, Function
    >>> f = Function('f')
    >>> y = f(x)
    >>> genform = a*y.diff(x) - (b*y**2 + c*y/x + d/x**2)
    >>> sol = dsolve(genform, y)
    >>> pprint(sol, wrap_line=False)
            /                                 /        __________________       \\
            |           __________________    |       /                2        ||
            |          /                2     |     \/  4*b*d - (a + c)  *log(x)||
           -|a + c - \/  4*b*d - (a + c)  *tan|C1 + ----------------------------||
            \                                 \                 2*a             //
    f(x) = ------------------------------------------------------------------------
                                            2*b*x

    >>> checkodesol(genform, sol, order=1)[0]
    True

    References
    ==========

    1. http://www.maplesoft.com/support/help/Maple/view.aspx?path=odeadvisor/Riccati
    2. http://eqworld.ipmnet.ru/en/solutions/ode/ode0106.pdf -
       http://eqworld.ipmnet.ru/en/solutions/ode/ode0123.pdf
    """
    hint = "Riccati_special_minus2"
    has_integral = False
    order = [1]

    def _wilds(self, f, x, order):
        a = Wild('a', exclude=[x, f(x), f(x).diff(x), 0])
        b = Wild('b', exclude=[x, f(x), f(x).diff(x), 0])
        c = Wild('c', exclude=[x, f(x), f(x).diff(x)])
        d = Wild('d', exclude=[x, f(x), f(x).diff(x)])
        return a, b, c, d

    def _equation(self, fx, x, order):
        a, b, c, d = self.wilds()
        return a*fx.diff(x) + b*fx**2 + c*fx/x + d/x**2

    def _get_general_solution(self, *, simplify_flag: bool = True):
        a, b, c, d = self.wilds_match()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        (C1,) = self.ode_problem.get_numbered_constants(num=1)
        mu = sqrt(4*d*b - (a - c)**2)

        gensol = Eq(fx, (a - c - mu*tan(mu/(2*a)*log(x) + C1))/(2*b*x))
        return [gensol]


class SecondNonlinearAutonomousConserved(SinglePatternODESolver):
    r"""
    Gives solution for the autonomous second order nonlinear
    differential equation of the form

    .. math :: f''(x) = g(f(x))

    The solution for this differential equation can be computed
    by multiplying by `f'(x)` and integrating on both sides,
    converting it into a first order differential equation.

    Examples
    ========

    >>> from sympy import Function, symbols, dsolve
    >>> f, g = symbols('f g', cls=Function)
    >>> x = symbols('x')

    >>> eq = f(x).diff(x, 2) - g(f(x))
    >>> dsolve(eq, simplify=False)
    [Eq(Integral(1/sqrt(C1 + 2*Integral(g(_u), _u)), (_u, f(x))), C2 + x),
    Eq(Integral(1/sqrt(C1 + 2*Integral(g(_u), _u)), (_u, f(x))), C2 - x)]

    >>> from sympy import exp, log
    >>> eq = f(x).diff(x, 2) - exp(f(x)) + log(f(x))
    >>> dsolve(eq, simplify=False)
    [Eq(Integral(1/sqrt(-2*_u*log(_u) + 2*_u + C1 + 2*exp(_u)), (_u, f(x))), C2 + x),
    Eq(Integral(1/sqrt(-2*_u*log(_u) + 2*_u + C1 + 2*exp(_u)), (_u, f(x))), C2 - x)]

    References
    ==========

    http://eqworld.ipmnet.ru/en/solutions/ode/ode0301.pdf
    """
    hint = "2nd_nonlinear_autonomous_conserved"
    has_integral = True
    order = [2]

    def _wilds(self, f, x, order):
        fy = Wild('fy', exclude=[0, f(x).diff(x), f(x).diff(x, 2)])
        return (fy,)

    def _equation(self, fx, x, order):
        fy = self.wilds()[0]
        return fx.diff(x, 2) + fy

    def _verify(self, fx):
        return self.ode_problem.is_autonomous

    def _get_general_solution(self, *, simplify_flag: bool = True):
        g = self.wilds_match()[0]
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        u = Dummy('u')
        g = g.subs(fx, u)
        C1, C2 = self.ode_problem.get_numbered_constants(num=2)
        inside = -2*Integral(g, u) + C1
        lhs = Integral(1/sqrt(inside), (u, fx))
        return [Eq(lhs, C2 + x), Eq(lhs, C2 - x)]


class Liouville(SinglePatternODESolver):
    r"""
    Solves 2nd order Liouville differential equations.

    The general form of a Liouville ODE is

    .. math:: \frac{d^2 y}{dx^2} + g(y) \left(\!
                \frac{dy}{dx}\!\right)^2 + h(x)
                \frac{dy}{dx}\text{.}

    The general solution is:

        >>> from sympy import Function, dsolve, Eq, pprint, diff
        >>> from sympy.abc import x
        >>> f, g, h = map(Function, ['f', 'g', 'h'])
        >>> genform = Eq(diff(f(x),x,x) + g(f(x))*diff(f(x),x)**2 +
        ... h(x)*diff(f(x),x), 0)
        >>> pprint(genform)
                          2                    2
                /d       \         d          d
        g(f(x))*|--(f(x))|  + h(x)*--(f(x)) + ---(f(x)) = 0
                \dx      /         dx           2
                                              dx
        >>> pprint(dsolve(genform, f(x), hint='Liouville_Integral'))
                                          f(x)
                  /                     /
                 |                     |
                 |     /               |     /
                 |    |                |    |
                 |  - | h(x) dx        |    | g(y) dy
                 |    |                |    |
                 |   /                 |   /
        C1 + C2* | e            dx +   |  e           dy = 0
                 |                     |
                /                     /

    Examples
    ========

    >>> from sympy import Function, dsolve, Eq, pprint
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(diff(f(x), x, x) + diff(f(x), x)**2/f(x) +
    ... diff(f(x), x)/x, f(x), hint='Liouville'))
               ________________           ________________
    [f(x) = -\/ C1 + C2*log(x) , f(x) = \/ C1 + C2*log(x) ]

    References
    ==========

    - Goldstein and Braun, "Advanced Methods for the Solution of Differential
      Equations", pp. 98
    - http://www.maplesoft.com/support/help/Maple/view.aspx?path=odeadvisor/Liouville

    # indirect doctest

    """
    hint = "Liouville"
    has_integral = True
    order = [2]

    def _wilds(self, f, x, order):
        d = Wild('d', exclude=[f(x).diff(x), f(x).diff(x, 2)])
        e = Wild('e', exclude=[f(x).diff(x)])
        k = Wild('k', exclude=[f(x).diff(x)])
        return d, e, k

    def _equation(self, fx, x, order):
        # Liouville ODE in the form
        # f(x).diff(x, 2) + g(f(x))*(f(x).diff(x))**2 + h(x)*f(x).diff(x)
        # See Goldstein and Braun, "Advanced Methods for the Solution of
        # Differential Equations", pg. 98
        d, e, k = self.wilds()
        return d*fx.diff(x, 2) + e*fx.diff(x)**2 + k*fx.diff(x)

    def _verify(self, fx):
        d, e, k = self.wilds_match()
        self.y = Dummy('y')
        x = self.ode_problem.sym
        self.g = simplify(e/d).subs(fx, self.y)
        self.h = simplify(k/d).subs(fx, self.y)
        if self.y in self.h.free_symbols or x in self.g.free_symbols:
            return False
        return True

    def _get_general_solution(self, *, simplify_flag: bool = True):
        d, e, k = self.wilds_match()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        C1, C2 = self.ode_problem.get_numbered_constants(num=2)
        int = Integral(exp(Integral(self.g, self.y)), (self.y, None, fx))
        gen_sol = Eq(int + C1*Integral(exp(-Integral(self.h, x)), x) + C2, 0)

        return [gen_sol]


class Separable(SinglePatternODESolver):
    r"""
    Solves separable 1st order differential equations.

    This is any differential equation that can be written as `P(y)
    \tfrac{dy}{dx} = Q(x)`.  The solution can then just be found by
    rearranging terms and integrating: `\int P(y) \,dy = \int Q(x) \,dx`.
    This hint uses :py:meth:`sympy.simplify.simplify.separatevars` as its back
    end, so if a separable equation is not caught by this solver, it is most
    likely the fault of that function.
    :py:meth:`~sympy.simplify.simplify.separatevars` is
    smart enough to do most expansion and factoring necessary to convert a
    separable equation `F(x, y)` into the proper form `P(x)\cdot{}Q(y)`.  The
    general solution is::

        >>> from sympy import Function, dsolve, Eq, pprint
        >>> from sympy.abc import x
        >>> a, b, c, d, f = map(Function, ['a', 'b', 'c', 'd', 'f'])
        >>> genform = Eq(a(x)*b(f(x))*f(x).diff(x), c(x)*d(f(x)))
        >>> pprint(genform)
                     d
        a(x)*b(f(x))*--(f(x)) = c(x)*d(f(x))
                     dx
        >>> pprint(dsolve(genform, f(x), hint='separable_Integral'))
             f(x)
           /                  /
          |                  |
          |  b(y)            | c(x)
          |  ---- dy = C1 +  | ---- dx
          |  d(y)            | a(x)
          |                  |
         /                  /

    Examples
    ========

    >>> from sympy import Function, dsolve, Eq
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(Eq(f(x)*f(x).diff(x) + x, 3*x*f(x)**2), f(x),
    ... hint='separable', simplify=False))
       /   2       \         2
    log\3*f (x) - 1/        x
    ---------------- = C1 + --
           6                2

    References
    ==========

    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 52

    # indirect doctest

    """
    hint = "separable"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        d = Wild('d', exclude=[f(x).diff(x), f(x).diff(x, 2)])
        e = Wild('e', exclude=[f(x).diff(x)])
        return d, e

    def _equation(self, fx, x, order):
        d, e = self.wilds()
        return d + e*fx.diff(x)

    def _verify(self, fx):
        d, e = self.wilds_match()
        self.y = Dummy('y')
        x = self.ode_problem.sym
        d = separatevars(d.subs(fx, self.y))
        e = separatevars(e.subs(fx, self.y))
        # m1[coeff]*m1[x]*m1[y] + m2[coeff]*m2[x]*m2[y]*y'
        self.m1 = separatevars(d, dict=True, symbols=(x, self.y))
        self.m2 = separatevars(e, dict=True, symbols=(x, self.y))
        if self.m1 and self.m2:
            return True
        return False

    def _get_match_object(self):
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        return self.m1, self.m2, x, fx

    def _get_general_solution(self, *, simplify_flag: bool = True):
        m1, m2, x, fx = self._get_match_object()
        (C1, ) = self.ode_problem.get_numbered_constants(num=1)
        int = Integral(m2['coeff']*m2[self.y]/m1[self.y],
        (self.y, None, fx))
        gen_sol = Eq(int, Integral(-m1['coeff']*m1[x]/
        m2[x], x) + C1)
        return [gen_sol]


class SeparableReduced(Separable):
    r"""
    Solves a differential equation that can be reduced to the separable form.

    The general form of this equation is

    .. math:: y' + (y/x) H(x^n y) = 0\text{}.

    This can be solved by substituting `u(y) = x^n y`.  The equation then
    reduces to the separable form `\frac{u'}{u (\mathrm{power} - H(u))} -
    \frac{1}{x} = 0`.

    The general solution is:

        >>> from sympy import Function, dsolve, pprint
        >>> from sympy.abc import x, n
        >>> f, g = map(Function, ['f', 'g'])
        >>> genform = f(x).diff(x) + (f(x)/x)*g(x**n*f(x))
        >>> pprint(genform)
                         / n     \
        d          f(x)*g\x *f(x)/
        --(f(x)) + ---------------
        dx                x
        >>> pprint(dsolve(genform, hint='separable_reduced'))
         n
        x *f(x)
          /
         |
         |         1
         |    ------------ dy = C1 + log(x)
         |    y*(n - g(y))
         |
         /

    See Also
    ========
    :obj:`sympy.solvers.ode.single.Separable`

    Examples
    ========

    >>> from sympy import Function, pprint
    >>> from sympy.solvers.ode.ode import dsolve
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> d = f(x).diff(x)
    >>> eq = (x - x**2*f(x))*d - f(x)
    >>> dsolve(eq, hint='separable_reduced')
    [Eq(f(x), (1 - sqrt(C1*x**2 + 1))/x), Eq(f(x), (sqrt(C1*x**2 + 1) + 1)/x)]
    >>> pprint(dsolve(eq, hint='separable_reduced'))
                   ___________            ___________
                  /     2                /     2
            1 - \/  C1*x  + 1          \/  C1*x  + 1  + 1
    [f(x) = ------------------, f(x) = ------------------]
                    x                          x

    References
    ==========

    - Joel Moses, "Symbolic Integration - The Stormy Decade", Communications
      of the ACM, Volume 14, Number 8, August 1971, pp. 558
    """
    hint = "separable_reduced"
    has_integral = True
    order = [1]

    def _degree(self, expr, x):
        # Made this function to calculate the degree of
        # x in an expression. If expr will be of form
        # x**p*y, (wheare p can be variables/rationals) then it
        # will return p.
        for val in expr:
            if val.has(x):
                if isinstance(val, Pow) and val.as_base_exp()[0] == x:
                    return (val.as_base_exp()[1])
                elif val == x:
                    return (val.as_base_exp()[1])
                else:
                    return self._degree(val.args, x)
        return 0

    def _powers(self, expr):
        # this function will return all the different relative power of x w.r.t f(x).
        # expr = x**p * f(x)**q then it will return {p/q}.
        pows = set()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        self.y = Dummy('y')
        if isinstance(expr, Add):
            exprs = expr.atoms(Add)
        elif isinstance(expr, Mul):
            exprs = expr.atoms(Mul)
        elif isinstance(expr, Pow):
            exprs = expr.atoms(Pow)
        else:
            exprs = {expr}

        for arg in exprs:
            if arg.has(x):
                _, u = arg.as_independent(x, fx)
                pow = self._degree((u.subs(fx, self.y), ), x)/self._degree((u.subs(fx, self.y), ), self.y)
                pows.add(pow)
        return pows

    def _verify(self, fx):
        num, den = self.wilds_match()
        x = self.ode_problem.sym
        factor = simplify(x/fx*num/den)
        # Try representing factor in terms of x^n*y
        # where n is lowest power of x in factor;
        # first remove terms like sqrt(2)*3 from factor.atoms(Mul)
        num, dem = factor.as_numer_denom()
        num = expand(num)
        dem = expand(dem)
        pows = self._powers(num)
        pows.update(self._powers(dem))
        pows = list(pows)
        if(len(pows)==1) and pows[0]!=zoo:
            self.t = Dummy('t')
            self.r2 = {'t': self.t}
            num = num.subs(x**pows[0]*fx, self.t)
            dem = dem.subs(x**pows[0]*fx, self.t)
            test = num/dem
            free = test.free_symbols
            if len(free) == 1 and free.pop() == self.t:
                self.r2.update({'power' : pows[0], 'u' : test})
                return True
            return False
        return False

    def _get_match_object(self):
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        u = self.r2['u'].subs(self.r2['t'], self.y)
        ycoeff = 1/(self.y*(self.r2['power'] - u))
        m1 = {self.y: 1, x: -1/x, 'coeff': 1}
        m2 = {self.y: ycoeff, x: 1, 'coeff': 1}
        return m1, m2, x, x**self.r2['power']*fx


class HomogeneousCoeffSubsDepDivIndep(SinglePatternODESolver):
    r"""
    Solves a 1st order differential equation with homogeneous coefficients
    using the substitution `u_1 = \frac{\text{<dependent
    variable>}}{\text{<independent variable>}}`.

    This is a differential equation

    .. math:: P(x, y) + Q(x, y) dy/dx = 0

    such that `P` and `Q` are homogeneous and of the same order.  A function
    `F(x, y)` is homogeneous of order `n` if `F(x t, y t) = t^n F(x, y)`.
    Equivalently, `F(x, y)` can be rewritten as `G(y/x)` or `H(x/y)`.  See
    also the docstring of :py:meth:`~sympy.solvers.ode.homogeneous_order`.

    If the coefficients `P` and `Q` in the differential equation above are
    homogeneous functions of the same order, then it can be shown that the
    substitution `y = u_1 x` (i.e. `u_1 = y/x`) will turn the differential
    equation into an equation separable in the variables `x` and `u`.  If
    `h(u_1)` is the function that results from making the substitution `u_1 =
    f(x)/x` on `P(x, f(x))` and `g(u_2)` is the function that results from the
    substitution on `Q(x, f(x))` in the differential equation `P(x, f(x)) +
    Q(x, f(x)) f'(x) = 0`, then the general solution is::

        >>> from sympy import Function, dsolve, pprint
        >>> from sympy.abc import x
        >>> f, g, h = map(Function, ['f', 'g', 'h'])
        >>> genform = g(f(x)/x) + h(f(x)/x)*f(x).diff(x)
        >>> pprint(genform)
         /f(x)\    /f(x)\ d
        g|----| + h|----|*--(f(x))
         \ x  /    \ x  / dx
        >>> pprint(dsolve(genform, f(x),
        ... hint='1st_homogeneous_coeff_subs_dep_div_indep_Integral'))
                       f(x)
                       ----
                        x
                         /
                        |
                        |       -h(u1)
        log(x) = C1 +   |  ---------------- d(u1)
                        |  u1*h(u1) + g(u1)
                        |
                       /

    Where `u_1 h(u_1) + g(u_1) \ne 0` and `x \ne 0`.

    See also the docstrings of
    :obj:`~sympy.solvers.ode.single.HomogeneousCoeffBest` and
    :obj:`~sympy.solvers.ode.single.HomogeneousCoeffSubsIndepDivDep`.

    Examples
    ========

    >>> from sympy import Function, dsolve
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x),
    ... hint='1st_homogeneous_coeff_subs_dep_div_indep', simplify=False))
                          /          3   \
                          |3*f(x)   f (x)|
                       log|------ + -----|
                          |  x         3 |
                          \           x  /
    log(x) = log(C1) - -------------------
                                3

    References
    ==========

    - https://en.wikipedia.org/wiki/Homogeneous_differential_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 59

    # indirect doctest

    """
    hint = "1st_homogeneous_coeff_subs_dep_div_indep"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        d = Wild('d', exclude=[f(x).diff(x), f(x).diff(x, 2)])
        e = Wild('e', exclude=[f(x).diff(x)])
        return d, e

    def _equation(self, fx, x, order):
        d, e = self.wilds()
        return d + e*fx.diff(x)

    def _verify(self, fx):
        self.d, self.e = self.wilds_match()
        self.y = Dummy('y')
        x = self.ode_problem.sym
        self.d = separatevars(self.d.subs(fx, self.y))
        self.e = separatevars(self.e.subs(fx, self.y))
        ordera = homogeneous_order(self.d, x, self.y)
        orderb = homogeneous_order(self.e, x, self.y)
        if ordera == orderb and ordera is not None:
            self.u = Dummy('u')
            if simplify((self.d + self.u*self.e).subs({x: 1, self.y: self.u})) != 0:
                return True
            return False
        return False

    def _get_match_object(self):
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        self.u1 = Dummy('u1')
        xarg = 0
        yarg = 0
        return [self.d, self.e, fx, x, self.u, self.u1, self.y, xarg, yarg]

    def _get_general_solution(self, *, simplify_flag: bool = True):
        d, e, fx, x, u, u1, y, xarg, yarg = self._get_match_object()
        (C1, ) = self.ode_problem.get_numbered_constants(num=1)
        int = Integral(
            (-e/(d + u1*e)).subs({x: 1, y: u1}),
            (u1, None, fx/x))
        sol = logcombine(Eq(log(x), int + log(C1)), force=True)
        gen_sol = sol.subs(fx, u).subs(((u, u - yarg), (x, x - xarg), (u, fx)))
        return [gen_sol]


class HomogeneousCoeffSubsIndepDivDep(SinglePatternODESolver):
    r"""
    Solves a 1st order differential equation with homogeneous coefficients
    using the substitution `u_2 = \frac{\text{<independent
    variable>}}{\text{<dependent variable>}}`.

    This is a differential equation

    .. math:: P(x, y) + Q(x, y) dy/dx = 0

    such that `P` and `Q` are homogeneous and of the same order.  A function
    `F(x, y)` is homogeneous of order `n` if `F(x t, y t) = t^n F(x, y)`.
    Equivalently, `F(x, y)` can be rewritten as `G(y/x)` or `H(x/y)`.  See
    also the docstring of :py:meth:`~sympy.solvers.ode.homogeneous_order`.

    If the coefficients `P` and `Q` in the differential equation above are
    homogeneous functions of the same order, then it can be shown that the
    substitution `x = u_2 y` (i.e. `u_2 = x/y`) will turn the differential
    equation into an equation separable in the variables `y` and `u_2`.  If
    `h(u_2)` is the function that results from making the substitution `u_2 =
    x/f(x)` on `P(x, f(x))` and `g(u_2)` is the function that results from the
    substitution on `Q(x, f(x))` in the differential equation `P(x, f(x)) +
    Q(x, f(x)) f'(x) = 0`, then the general solution is:

    >>> from sympy import Function, dsolve, pprint
    >>> from sympy.abc import x
    >>> f, g, h = map(Function, ['f', 'g', 'h'])
    >>> genform = g(x/f(x)) + h(x/f(x))*f(x).diff(x)
    >>> pprint(genform)
     / x  \    / x  \ d
    g|----| + h|----|*--(f(x))
     \f(x)/    \f(x)/ dx
    >>> pprint(dsolve(genform, f(x),
    ... hint='1st_homogeneous_coeff_subs_indep_div_dep_Integral'))
                 x
                ----
                f(x)
                  /
                 |
                 |       -g(u1)
                 |  ---------------- d(u1)
                 |  u1*g(u1) + h(u1)
                 |
                /
    <BLANKLINE>
    f(x) = C1*e

    Where `u_1 g(u_1) + h(u_1) \ne 0` and `f(x) \ne 0`.

    See also the docstrings of
    :obj:`~sympy.solvers.ode.single.HomogeneousCoeffBest` and
    :obj:`~sympy.solvers.ode.single.HomogeneousCoeffSubsDepDivIndep`.

    Examples
    ========

    >>> from sympy import Function, pprint, dsolve
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x),
    ... hint='1st_homogeneous_coeff_subs_indep_div_dep',
    ... simplify=False))
                             /    2    \
                             | 3*x     |
                          log|----- + 1|
                             | 2       |
                             \f (x)    /
    log(f(x)) = log(C1) - --------------
                                3

    References
    ==========

    - https://en.wikipedia.org/wiki/Homogeneous_differential_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 59

    # indirect doctest

    """
    hint = "1st_homogeneous_coeff_subs_indep_div_dep"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        d = Wild('d', exclude=[f(x).diff(x), f(x).diff(x, 2)])
        e = Wild('e', exclude=[f(x).diff(x)])
        return d, e

    def _equation(self, fx, x, order):
        d, e = self.wilds()
        return d + e*fx.diff(x)

    def _verify(self, fx):
        self.d, self.e = self.wilds_match()
        self.y = Dummy('y')
        x = self.ode_problem.sym
        self.d = separatevars(self.d.subs(fx, self.y))
        self.e = separatevars(self.e.subs(fx, self.y))
        ordera = homogeneous_order(self.d, x, self.y)
        orderb = homogeneous_order(self.e, x, self.y)
        if ordera == orderb and ordera is not None:
            self.u = Dummy('u')
            if simplify((self.e + self.u*self.d).subs({x: self.u, self.y: 1})) != 0:
                return True
            return False
        return False

    def _get_match_object(self):
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        self.u1 = Dummy('u1')
        xarg = 0
        yarg = 0
        return [self.d, self.e, fx, x, self.u, self.u1, self.y, xarg, yarg]

    def _get_general_solution(self, *, simplify_flag: bool = True):
        d, e, fx, x, u, u1, y, xarg, yarg = self._get_match_object()
        (C1, ) = self.ode_problem.get_numbered_constants(num=1)
        int = Integral(simplify((-d/(e + u1*d)).subs({x: u1, y: 1})),(u1, None, x/fx))
        sol = logcombine(Eq(log(fx), int + log(C1)), force=True)
        gen_sol = sol.subs(fx, u).subs(((u, u - yarg), (x, x - xarg), (u, fx)))
        return [gen_sol]


class HomogeneousCoeffBest(HomogeneousCoeffSubsIndepDivDep, HomogeneousCoeffSubsDepDivIndep):
    r"""
    Returns the best solution to an ODE from the two hints
    ``1st_homogeneous_coeff_subs_dep_div_indep`` and
    ``1st_homogeneous_coeff_subs_indep_div_dep``.

    This is as determined by :py:meth:`~sympy.solvers.ode.ode.ode_sol_simplicity`.

    See the
    :obj:`~sympy.solvers.ode.single.HomogeneousCoeffSubsIndepDivDep`
    and
    :obj:`~sympy.solvers.ode.single.HomogeneousCoeffSubsDepDivIndep`
    docstrings for more information on these hints.  Note that there is no
    ``ode_1st_homogeneous_coeff_best_Integral`` hint.

    Examples
    ========

    >>> from sympy import Function, dsolve, pprint
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x),
    ... hint='1st_homogeneous_coeff_best', simplify=False))
                             /    2    \
                             | 3*x     |
                          log|----- + 1|
                             | 2       |
                             \f (x)    /
    log(f(x)) = log(C1) - --------------
                                3

    References
    ==========

    - https://en.wikipedia.org/wiki/Homogeneous_differential_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 59

    # indirect doctest

    """
    hint = "1st_homogeneous_coeff_best"
    has_integral = False
    order = [1]

    def _verify(self, fx):
        if HomogeneousCoeffSubsIndepDivDep._verify(self, fx) and HomogeneousCoeffSubsDepDivIndep._verify(self, fx):
            return True
        return False

    def _get_general_solution(self, *, simplify_flag: bool = True):
        # There are two substitutions that solve the equation, u1=y/x and u2=x/y
        # # They produce different integrals, so try them both and see which
        # # one is easier
        sol1 = HomogeneousCoeffSubsIndepDivDep._get_general_solution(self)
        sol2 = HomogeneousCoeffSubsDepDivIndep._get_general_solution(self)
        fx = self.ode_problem.func
        if simplify_flag:
            sol1 = odesimp(self.ode_problem.eq, *sol1, fx, "1st_homogeneous_coeff_subs_indep_div_dep")
            sol2 = odesimp(self.ode_problem.eq, *sol2, fx, "1st_homogeneous_coeff_subs_dep_div_indep")
        return min([sol1, sol2], key=lambda x: ode_sol_simplicity(x, fx, trysolving=not simplify))


class LinearCoefficients(HomogeneousCoeffBest):
    r"""
    Solves a differential equation with linear coefficients.

    The general form of a differential equation with linear coefficients is

    .. math:: y' + F\left(\!\frac{a_1 x + b_1 y + c_1}{a_2 x + b_2 y +
                c_2}\!\right) = 0\text{,}

    where `a_1`, `b_1`, `c_1`, `a_2`, `b_2`, `c_2` are constants and `a_1 b_2
    - a_2 b_1 \ne 0`.

    This can be solved by substituting:

    .. math:: x = x' + \frac{b_2 c_1 - b_1 c_2}{a_2 b_1 - a_1 b_2}

              y = y' + \frac{a_1 c_2 - a_2 c_1}{a_2 b_1 - a_1
                  b_2}\text{.}

    This substitution reduces the equation to a homogeneous differential
    equation.

    See Also
    ========
    :obj:`sympy.solvers.ode.single.HomogeneousCoeffBest`
    :obj:`sympy.solvers.ode.single.HomogeneousCoeffSubsIndepDivDep`
    :obj:`sympy.solvers.ode.single.HomogeneousCoeffSubsDepDivIndep`

    Examples
    ========

    >>> from sympy import Function, pprint
    >>> from sympy.solvers.ode.ode import dsolve
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> df = f(x).diff(x)
    >>> eq = (x + f(x) + 1)*df + (f(x) - 6*x + 1)
    >>> dsolve(eq, hint='linear_coefficients')
    [Eq(f(x), -x - sqrt(C1 + 7*x**2) - 1), Eq(f(x), -x + sqrt(C1 + 7*x**2) - 1)]
    >>> pprint(dsolve(eq, hint='linear_coefficients'))
                      ___________                     ___________
                   /         2                     /         2
    [f(x) = -x - \/  C1 + 7*x   - 1, f(x) = -x + \/  C1 + 7*x   - 1]


    References
    ==========

    - Joel Moses, "Symbolic Integration - The Stormy Decade", Communications
      of the ACM, Volume 14, Number 8, August 1971, pp. 558
    """
    hint = "linear_coefficients"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        d = Wild('d', exclude=[f(x).diff(x), f(x).diff(x, 2)])
        e = Wild('e', exclude=[f(x).diff(x)])
        return d, e

    def _equation(self, fx, x, order):
        d, e = self.wilds()
        return d + e*fx.diff(x)

    def _verify(self, fx):
        self.d, self.e = self.wilds_match()
        a, b = self.wilds()
        F = self.d/self.e
        x = self.ode_problem.sym
        params = self._linear_coeff_match(F, fx)
        if params:
            self.xarg, self.yarg = params
            u = Dummy('u')
            t = Dummy('t')
            self.y = Dummy('y')
            # Dummy substitution for df and f(x).
            dummy_eq = self.ode_problem.eq.subs(((fx.diff(x), t), (fx, u)))
            reps = ((x, x + self.xarg), (u, u + self.yarg), (t, fx.diff(x)), (u, fx))
            dummy_eq = simplify(dummy_eq.subs(reps))
            # get the re-cast values for e and d
            r2 = collect(expand(dummy_eq), [fx.diff(x), fx]).match(a*fx.diff(x) + b)
            if r2:
                self.d, self.e = r2[b], r2[a]
                orderd = homogeneous_order(self.d, x, fx)
                ordere = homogeneous_order(self.e, x, fx)
                if orderd == ordere and orderd is not None:
                    self.d = self.d.subs(fx, self.y)
                    self.e = self.e.subs(fx, self.y)
                    return True
                return False
            return False

    def _linear_coeff_match(self,expr, func):
        r"""
        Helper function to match hint ``linear_coefficients``.

        Matches the expression to the form `(a_1 x + b_1 f(x) + c_1)/(a_2 x + b_2
        f(x) + c_2)` where the following conditions hold:

        1. `a_1`, `b_1`, `c_1`, `a_2`, `b_2`, `c_2` are Rationals;
        2. `c_1` or `c_2` are not equal to zero;
        3. `a_2 b_1 - a_1 b_2` is not equal to zero.

        Return ``xarg``, ``yarg`` where

        1. ``xarg`` = `(b_2 c_1 - b_1 c_2)/(a_2 b_1 - a_1 b_2)`
        2. ``yarg`` = `(a_1 c_2 - a_2 c_1)/(a_2 b_1 - a_1 b_2)`


        Examples
        ========

        >>> from sympy import Function
        >>> from sympy.abc import x
        >>> from sympy.solvers.ode.single import LinearCoefficients
        >>> from sympy.functions.elementary.trigonometric import sin
        >>> f = Function('f')
        >>> eq = (-25*f(x) - 8*x + 62)/(4*f(x) + 11*x - 11)
        >>> obj = LinearCoefficients(eq)
        >>> obj._linear_coeff_match(eq, f(x))
        (1/9, 22/9)
        >>> eq = sin((-5*f(x) - 8*x + 6)/(4*f(x) + x - 1))
        >>> obj = LinearCoefficients(eq)
        >>> obj._linear_coeff_match(eq, f(x))
        (19/27, 2/27)
        >>> eq = sin(f(x)/x)
        >>> obj = LinearCoefficients(eq)
        >>> obj._linear_coeff_match(eq, f(x))

        """
        f = func.func
        x = func.args[0]
        def abc(eq):
            r'''
            Internal function of _linear_coeff_match
            that returns Rationals a, b, c
            if eq is a*x + b*f(x) + c, else None.
            '''
            eq = _mexpand(eq)
            c = eq.as_independent(x, f(x), as_Add=True)[0]
            if not c.is_Rational:
                return
            a = eq.coeff(x)
            if not a.is_Rational:
                return
            b = eq.coeff(f(x))
            if not b.is_Rational:
                return
            if eq == a*x + b*f(x) + c:
                return a, b, c

        def match(arg):
            r'''
            Internal function of _linear_coeff_match that returns Rationals a1,
            b1, c1, a2, b2, c2 and a2*b1 - a1*b2 of the expression (a1*x + b1*f(x)
            + c1)/(a2*x + b2*f(x) + c2) if one of c1 or c2 and a2*b1 - a1*b2 is
            non-zero, else None.
            '''
            n, d = arg.together().as_numer_denom()
            m = abc(n)
            if m is not None:
                a1, b1, c1 = m
                m = abc(d)
                if m is not None:
                    a2, b2, c2 = m
                    d = a2*b1 - a1*b2
                    if (c1 or c2) and d:
                        return a1, b1, c1, a2, b2, c2, d

        m = [fi.args[0] for fi in expr.atoms(Function) if fi.func != f and
            len(fi.args) == 1 and not fi.args[0].is_Function] or {expr}
        m1 = match(m.pop())
        if m1 and all(match(mi) == m1 for mi in m):
            a1, b1, c1, a2, b2, c2, denom = m1
            return (b2*c1 - b1*c2)/denom, (a1*c2 - a2*c1)/denom

    def _get_match_object(self):
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        self.u1 = Dummy('u1')
        u = Dummy('u')
        return [self.d, self.e, fx, x, u, self.u1, self.y, self.xarg, self.yarg]


class NthOrderReducible(SingleODESolver):
    r"""
    Solves ODEs that only involve derivatives of the dependent variable using
    a substitution of the form `f^n(x) = g(x)`.

    For example any second order ODE of the form `f''(x) = h(f'(x), x)` can be
    transformed into a pair of 1st order ODEs `g'(x) = h(g(x), x)` and
    `f'(x) = g(x)`. Usually the 1st order ODE for `g` is easier to solve. If
    that gives an explicit solution for `g` then `f` is found simply by
    integration.


    Examples
    ========

    >>> from sympy import Function, dsolve, Eq
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> eq = Eq(x*f(x).diff(x)**2 + f(x).diff(x, 2), 0)
    >>> dsolve(eq, f(x), hint='nth_order_reducible')
    ... # doctest: +NORMALIZE_WHITESPACE
    Eq(f(x), C1 - sqrt(-1/C2)*log(-C2*sqrt(-1/C2) + x) + sqrt(-1/C2)*log(C2*sqrt(-1/C2) + x))

    """
    hint = "nth_order_reducible"
    has_integral = False

    def _matches(self):
        # Any ODE that can be solved with a substitution and
        # repeated integration e.g.:
        # `d^2/dx^2(y) + x*d/dx(y) = constant
        #f'(x) must be finite for this to work
        eq = self.ode_problem.eq_preprocessed
        func = self.ode_problem.func
        x = self.ode_problem.sym
        r"""
        Matches any differential equation that can be rewritten with a smaller
        order. Only derivatives of ``func`` alone, wrt a single variable,
        are considered, and only in them should ``func`` appear.
        """
        # ODE only handles functions of 1 variable so this affirms that state
        assert len(func.args) == 1
        vc = [d.variable_count[0] for d in eq.atoms(Derivative)
            if d.expr == func and len(d.variable_count) == 1]
        ords = [c for v, c in vc if v == x]
        if len(ords) < 2:
            return False
        self.smallest = min(ords)
        # make sure func does not appear outside of derivatives
        D = Dummy()
        if eq.subs(func.diff(x, self.smallest), D).has(func):
            return False
        return True

    def _get_general_solution(self, *, simplify_flag: bool = True):
        eq = self.ode_problem.eq
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        n = self.smallest
        # get a unique function name for g
        names = [a.name for a in eq.atoms(AppliedUndef)]
        while True:
            name = Dummy().name
            if name not in names:
                g = Function(name)
                break
        w = f(x).diff(x, n)
        geq = eq.subs(w, g(x))
        gsol = dsolve(geq, g(x))

        if not isinstance(gsol, list):
            gsol = [gsol]

        # Might be multiple solutions to the reduced ODE:
        fsol = []
        for gsoli in gsol:
            fsoli = dsolve(gsoli.subs(g(x), w), f(x))  # or do integration n times
            fsol.append(fsoli)

        return fsol


class Hypergeometric2nd(SingleODESolver):
    r"""
    Solves 2nd order linear differential equations.

    It computes special function solutions which can be expressed using the
    2F1, 1F1 or 0F1 hypergeometric functions.

    .. math:: y'' + A(x) y' + B(x) y = 0\text{,}

    where `A` and `B` are rational functions.

    These kinds of differential equations have solution of non-Liouvillian form.

    Given linear ODE can be obtained from 2F1 given by

    .. math:: (x^2 - x) y'' + ((a + b + 1) x - c) y' + b a y = 0\text{,}

    where {a, b, c} are arbitrary constants.

    Notes
    =====

    The algorithm should find any solution of the form

    .. math:: y = P(x) _pF_q(..; ..;\frac{\alpha x^k + \beta}{\gamma x^k + \delta})\text{,}

    where pFq is any of 2F1, 1F1 or 0F1 and `P` is an "arbitrary function".
    Currently only the 2F1 case is implemented in SymPy but the other cases are
    described in the paper and could be implemented in future (contributions
    welcome!).


    Examples
    ========

    >>> from sympy import Function, dsolve, pprint
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> eq = (x*x - x)*f(x).diff(x,2) + (5*x - 1)*f(x).diff(x) + 4*f(x)
    >>> pprint(dsolve(eq, f(x), '2nd_hypergeometric'))
                                        _
           /        /           4  \\  |_  /-1, -1 |  \
           |C1 + C2*|log(x) + -----||* |   |       | x|
           \        \         x + 1// 2  1 \  1    |  /
    f(x) = --------------------------------------------
                                    3
                             (x - 1)


    References
    ==========

    - "Non-Liouvillian solutions for second order linear ODEs" by L. Chan, E.S. Cheb-Terrab

    """
    hint = "2nd_hypergeometric"
    has_integral = True

    def _matches(self):
        eq = self.ode_problem.eq_preprocessed
        func = self.ode_problem.func
        r = match_2nd_hypergeometric(eq, func)
        self.match_object = None
        if r:
            A, B = r
            d = equivalence_hypergeometric(A, B, func)
            if d:
                if d['type'] == "2F1":
                    self.match_object = match_2nd_2F1_hypergeometric(d['I0'], d['k'], d['sing_point'], func)
                    if self.match_object is not None:
                        self.match_object.update({'A':A, 'B':B})
            # We can extend it for 1F1 and 0F1 type also.
        return self.match_object is not None

    def _get_general_solution(self, *, simplify_flag: bool = True):
        eq = self.ode_problem.eq
        func = self.ode_problem.func
        if self.match_object['type'] == "2F1":
            sol = get_sol_2F1_hypergeometric(eq, func, self.match_object)
            if sol is None:
                raise NotImplementedError("The given ODE " + str(eq) + " cannot be solved by"
                    + " the hypergeometric method")

        return [sol]


class NthConstCoeffHomogen(SingleODESolver):
    r"""
    Solves an `n`\th order linear homogeneous differential equation with
    constant coefficients.

    This is an equation of the form

    .. math:: a_n f^{(n)}(x) + a_{n-1} f^{(n-1)}(x) + \cdots + a_1 f'(x)
                + a_0 f(x) = 0\text{.}

    These equations can be solved in a general manner, by taking the roots of
    the characteristic equation `a_n m^n + a_{n-1} m^{n-1} + \cdots + a_1 m +
    a_0 = 0`.  The solution will then be the sum of `C_n x^i e^{r x}` terms,
    for each where `C_n` is an arbitrary constant, `r` is a root of the
    characteristic equation and `i` is one of each from 0 to the multiplicity
    of the root - 1 (for example, a root 3 of multiplicity 2 would create the
    terms `C_1 e^{3 x} + C_2 x e^{3 x}`).  The exponential is usually expanded
    for complex roots using Euler's equation `e^{I x} = \cos(x) + I \sin(x)`.
    Complex roots always come in conjugate pairs in polynomials with real
    coefficients, so the two roots will be represented (after simplifying the
    constants) as `e^{a x} \left(C_1 \cos(b x) + C_2 \sin(b x)\right)`.

    If SymPy cannot find exact roots to the characteristic equation, a
    :py:class:`~sympy.polys.rootoftools.ComplexRootOf` instance will be return
    instead.

    >>> from sympy import Function, dsolve
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> dsolve(f(x).diff(x, 5) + 10*f(x).diff(x) - 2*f(x), f(x),
    ... hint='nth_linear_constant_coeff_homogeneous')
    ... # doctest: +NORMALIZE_WHITESPACE
    Eq(f(x), C5*exp(x*CRootOf(_x**5 + 10*_x - 2, 0))
    + (C1*sin(x*im(CRootOf(_x**5 + 10*_x - 2, 1)))
    + C2*cos(x*im(CRootOf(_x**5 + 10*_x - 2, 1))))*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 1)))
    + (C3*sin(x*im(CRootOf(_x**5 + 10*_x - 2, 3)))
    + C4*cos(x*im(CRootOf(_x**5 + 10*_x - 2, 3))))*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 3))))

    Note that because this method does not involve integration, there is no
    ``nth_linear_constant_coeff_homogeneous_Integral`` hint.

    Examples
    ========

    >>> from sympy import Function, dsolve, pprint
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(f(x).diff(x, 4) + 2*f(x).diff(x, 3) -
    ... 2*f(x).diff(x, 2) - 6*f(x).diff(x) + 5*f(x), f(x),
    ... hint='nth_linear_constant_coeff_homogeneous'))
                        x                            -2*x
    f(x) = (C1 + C2*x)*e  + (C3*sin(x) + C4*cos(x))*e

    References
    ==========

    - https://en.wikipedia.org/wiki/Linear_differential_equation section:
      Nonhomogeneous_equation_with_constant_coefficients
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 211

    # indirect doctest

    """
    hint = "nth_linear_constant_coeff_homogeneous"
    has_integral = False

    def _matches(self):
        eq = self.ode_problem.eq_high_order_free
        func = self.ode_problem.func
        order = self.ode_problem.order
        x = self.ode_problem.sym
        self.r = self.ode_problem.get_linear_coefficients(eq, func, order)
        if order and self.r and not any(self.r[i].has(x) for i in self.r if i >= 0):
            if not self.r[-1]:
                return True
            else:
                return False
        return False

    def _get_sols(self,r):
        x = self.ode_problem.sym
        order = self.ode_problem.order
        # First, set up characteristic equation.
        chareq, symbol = S.Zero, Dummy('x')

        for i in r.keys():
            if type(i) == str or i < 0:
                pass
            else:
                chareq += r[i]*symbol**i

        chareq = Poly(chareq, symbol)
        # Can't just call roots because it doesn't return rootof for unsolveable
        # polynomials.
        chareqroots = roots(chareq, multiple=True)
        if len(chareqroots) != order:
            chareqroots = [rootof(chareq, k) for k in range(chareq.degree())]

        chareq_is_complex = not all([i.is_real for i in chareq.all_coeffs()])

        # Create a dict root: multiplicity or charroots
        charroots = defaultdict(int)
        for root in chareqroots:
            charroots[root] += 1
        # We need to keep track of terms so we can run collect() at the end.
        # This is necessary for constantsimp to work properly.
        collectterms = []
        gensols = []
        conjugate_roots = [] # used to prevent double-use of conjugate roots
        # Loop over roots in theorder provided by roots/rootof...
        for root in chareqroots:
            # but don't repoeat multiple roots.
            if root not in charroots:
                continue
            multiplicity = charroots.pop(root)
            for i in range(multiplicity):
                if chareq_is_complex:
                    gensols.append(x**i*exp(root*x))
                    collectterms = [(i, root, 0)] + collectterms
                    continue
                reroot = re(root)
                imroot = im(root)
                if imroot.has(atan2) and reroot.has(atan2):
                    # Remove this condition when re and im stop returning
                    # circular atan2 usages.
                    gensols.append(x**i*exp(root*x))
                    collectterms = [(i, root, 0)] + collectterms
                else:
                    if root in conjugate_roots:
                        collectterms = [(i, reroot, imroot)] + collectterms
                        continue
                    if imroot == 0:
                        gensols.append(x**i*exp(reroot*x))
                        collectterms = [(i, reroot, 0)] + collectterms
                        continue
                    conjugate_roots.append(conjugate(root))
                    gensols.append(x**i*exp(reroot*x) * sin(abs(imroot) * x))
                    gensols.append(x**i*exp(reroot*x) * cos(    imroot  * x))

                    # This ordering is important
                    collectterms = [(i, reroot, imroot)] + collectterms
        return gensols, collectterms

    def _get_simplified_sol(self, sol, collectterms):
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        collectterms.sort(key=default_sort_key)
        collectterms.reverse()
        assert len(sol) == 1 and sol[0].lhs == f(x)
        sol = sol[0].rhs
        sol = expand_mul(sol)
        for i, reroot, imroot in collectterms:
            sol = collect(sol, x**i*exp(reroot*x)*sin(abs(imroot)*x))
            sol = collect(sol, x**i*exp(reroot*x)*cos(imroot*x))
        for i, reroot, imroot in collectterms:
            sol = collect(sol, x**i*exp(reroot*x))
        sol = powsimp(sol)
        return [Eq(f(x), sol)]

    def _get_general_solution(self, *, simplify_flag: bool = True):
        fx = self.ode_problem.func
        gensols, collectterms = self._get_sols(self.r)
        # A generator of constants
        constants = self.ode_problem.get_numbered_constants(num=len(gensols))
        gsol = Add(*[i*j for (i, j) in zip(constants, gensols)])
        gsol = [Eq(fx, gsol)]
        if simplify_flag:
            gsol = self._get_simplified_sol(gsol, collectterms)

        return gsol


class NthConstCoeffVar(SingleODESolver):
    r"""
    Solves an `n`\th order linear differential equation with constant
    coefficients using the method of variation of parameters.

    This method works on any differential equations of the form

    .. math:: f^{(n)}(x) + a_{n-1} f^{(n-1)}(x) + \cdots + a_1 f'(x) + a_0
                f(x) = P(x)\text{.}

    This method works by assuming that the particular solution takes the form

    .. math:: \sum_{x=1}^{n} c_i(x) y_i(x)\text{,}

    where `y_i` is the `i`\th solution to the homogeneous equation.  The
    solution is then solved using Wronskian's and Cramer's Rule.  The
    particular solution is given by

    .. math:: \sum_{x=1}^n \left( \int \frac{W_i(x)}{W(x)} \,dx
                \right) y_i(x) \text{,}

    where `W(x)` is the Wronskian of the fundamental system (the system of `n`
    linearly independent solutions to the homogeneous equation), and `W_i(x)`
    is the Wronskian of the fundamental system with the `i`\th column replaced
    with `[0, 0, \cdots, 0, P(x)]`.

    This method is general enough to solve any `n`\th order inhomogeneous
    linear differential equation with constant coefficients, but sometimes
    SymPy cannot simplify the Wronskian well enough to integrate it.  If this
    method hangs, try using the
    ``nth_linear_constant_coeff_variation_of_parameters_Integral`` hint and
    simplifying the integrals manually.  Also, prefer using
    ``nth_linear_constant_coeff_undetermined_coefficients`` when it
    applies, because it doesn't use integration, making it faster and more
    reliable.

    Warning, using simplify=False with
    'nth_linear_constant_coeff_variation_of_parameters' in
    :py:meth:`~sympy.solvers.ode.dsolve` may cause it to hang, because it will
    not attempt to simplify the Wronskian before integrating.  It is
    recommended that you only use simplify=False with
    'nth_linear_constant_coeff_variation_of_parameters_Integral' for this
    method, especially if the solution to the homogeneous equation has
    trigonometric functions in it.

    Examples
    ========

    >>> from sympy import Function, dsolve, pprint, exp, log
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(f(x).diff(x, 3) - 3*f(x).diff(x, 2) +
    ... 3*f(x).diff(x) - f(x) - exp(x)*log(x), f(x),
    ... hint='nth_linear_constant_coeff_variation_of_parameters'))
           /       /       /     x*log(x)   11*x\\\  x
    f(x) = |C1 + x*|C2 + x*|C3 + -------- - ----|||*e
           \       \       \        6        36 ///

    References
    ==========

    - https://en.wikipedia.org/wiki/Variation_of_parameters
    - http://planetmath.org/VariationOfParameters
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 233

    # indirect doctest

    """
    hint = "nth_linear_constant_coeff_variation_of_parameters"
    has_integral = True

    def _matches(self):
        eq = self.ode_problem.eq_high_order_free
        func = self.ode_problem.func
        order = self.ode_problem.order
        x = self.ode_problem.sym
        self.r = self.ode_problem.get_linear_coefficients(eq, func, order)

        if order and self.r and not any(self.r[i].has(x) for i in self.r if i >= 0):
            if self.r[-1]:
                return True
            else:
                return False
        return False

    def _solve_variation_of_parameters(self,match_object):
        r"""
        Helper function for the method of variation of parameters and nonhomogeneous euler eq.

        See the
        :py:meth:`~sympy.solvers.ode.single.NthConstCoeffVar`
        docstring for more information on this method.

        The parameter ``match`` should be a dictionary that has the following
        keys:

        ``list``
        A list of solutions to the homogeneous equation.

        ``sol``
        The general solution.

        """
        eq = self.ode_problem.eq_high_order_free
        f = self.ode_problem.func.func
        order = self.ode_problem.order
        x = self.ode_problem.sym
        r = match_object
        psol = 0
        gensols = r['list']
        gsol = r['sol']
        wr = wronskian(gensols, x)

        if r.get('simplify_flag', True):
            wr = simplify(wr)  # We need much better simplification for
                            # some ODEs. See issue 4662, for example.
            # To reduce commonly occurring sin(x)**2 + cos(x)**2 to 1
            wr = trigsimp(wr, deep=True, recursive=True)
        if not wr:
            # The wronskian will be 0 iff the solutions are not linearly
            # independent.
            raise NotImplementedError("Cannot find " + str(order) +
            " solutions to the homogeneous equation necessary to apply " +
            "variation of parameters to " + str(eq) + " (Wronskian == 0)")
        if len(gensols) != order:
            raise NotImplementedError("Cannot find " + str(order) +
            " solutions to the homogeneous equation necessary to apply " +
            "variation of parameters to " +
            str(eq) + " (number of terms != order)")
        negoneterm = (-1)**(order)
        for i in gensols:
            psol += negoneterm*Integral(wronskian([sol for sol in gensols if sol != i], x)*r[-1]/wr, x)*i/r[order]
            negoneterm *= -1

        if r.get('simplify_flag', True):
            psol = simplify(psol)
            psol = trigsimp(psol, deep=True)
        return Eq(f(x), gsol.rhs + psol)

    def _get_general_solution(self, *, simplify_flag: bool = True):
        eq = self.ode_problem.eq_high_order_free
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        homogen_instance = NthConstCoeffHomogen(SingleODEProblem(eq, f(x), x))
        gensols, collectterms = homogen_instance._get_sols(self.r)
        # A generator of constants
        constants = self.ode_problem.get_numbered_constants(num=len(gensols))
        gsol = Add(*[i*j for (i, j) in zip(constants, gensols)])
        gsol = Eq(f(x), gsol)
        self.r.update({'list': gensols, 'sol': gsol, 'simpliy_flag': simplify_flag})
        gsol = self._solve_variation_of_parameters(self.r)
        if simplify_flag:
            gsol = homogen_instance._get_simplified_sol([gsol], collectterms)
        return gsol


class NthConstCoeffUndet(SingleODESolver):
    r"""
    Solves an `n`\th order linear differential equation with constant
    coefficients using the method of undetermined coefficients.

    This method works on differential equations of the form

    .. math:: a_n f^{(n)}(x) + a_{n-1} f^{(n-1)}(x) + \cdots + a_1 f'(x)
                + a_0 f(x) = P(x)\text{,}

    where `P(x)` is a function that has a finite number of linearly
    independent derivatives.

    Functions that fit this requirement are finite sums functions of the form
    `a x^i e^{b x} \sin(c x + d)` or `a x^i e^{b x} \cos(c x + d)`, where `i`
    is a non-negative integer and `a`, `b`, `c`, and `d` are constants.  For
    example any polynomial in `x`, functions like `x^2 e^{2 x}`, `x \sin(x)`,
    and `e^x \cos(x)` can all be used.  Products of `\sin`'s and `\cos`'s have
    a finite number of derivatives, because they can be expanded into `\sin(a
    x)` and `\cos(b x)` terms.  However, SymPy currently cannot do that
    expansion, so you will need to manually rewrite the expression in terms of
    the above to use this method.  So, for example, you will need to manually
    convert `\sin^2(x)` into `(1 + \cos(2 x))/2` to properly apply the method
    of undetermined coefficients on it.

    This method works by creating a trial function from the expression and all
    of its linear independent derivatives and substituting them into the
    original ODE.  The coefficients for each term will be a system of linear
    equations, which are be solved for and substituted, giving the solution.
    If any of the trial functions are linearly dependent on the solution to
    the homogeneous equation, they are multiplied by sufficient `x` to make
    them linearly independent.

    Examples
    ========

    >>> from sympy import Function, dsolve, pprint, exp, cos
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> pprint(dsolve(f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) -
    ... 4*exp(-x)*x**2 + cos(2*x), f(x),
    ... hint='nth_linear_constant_coeff_undetermined_coefficients'))
           /       /      3\\
           |       |     x ||  -x   4*sin(2*x)   3*cos(2*x)
    f(x) = |C1 + x*|C2 + --||*e   - ---------- + ----------
           \       \     3 //           25           25

    References
    ==========

    - https://en.wikipedia.org/wiki/Method_of_undetermined_coefficients
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 221

    # indirect doctest

    """
    hint = "nth_linear_constant_coeff_undetermined_coefficients"
    has_integral = False

    def _matches(self):
        eq = self.ode_problem.eq_high_order_free
        func = self.ode_problem.func
        order = self.ode_problem.order
        x = self.ode_problem.sym
        self.r = self.ode_problem.get_linear_coefficients(eq, func, order)
        does_match = False
        if order and self.r and not any(self.r[i].has(x) for i in self.r if i >= 0):
            if self.r[-1]:
                eq_homogeneous = Add(eq,-self.r[-1])
                undetcoeff = self._undetermined_coefficients_match(self.r[-1], x, func, eq_homogeneous)
                if undetcoeff['test']:
                    self.r['trialset'] = undetcoeff['trialset']
                    does_match = True
        return does_match

    def _undetermined_coefficients_match(self,expr, x, func=None, eq_homogeneous=S.Zero):
        r"""
        Returns a trial function match if undetermined coefficients can be applied
        to ``expr``, and ``None`` otherwise.

        A trial expression can be found for an expression for use with the method
        of undetermined coefficients if the expression is an
        additive/multiplicative combination of constants, polynomials in `x` (the
        independent variable of expr), `\sin(a x + b)`, `\cos(a x + b)`, and
        `e^{a x}` terms (in other words, it has a finite number of linearly
        independent derivatives).

        Note that you may still need to multiply each term returned here by
        sufficient `x` to make it linearly independent with the solutions to the
        homogeneous equation.

        This is intended for internal use by ``undetermined_coefficients`` hints.

        SymPy currently has no way to convert `\sin^n(x) \cos^m(y)` into a sum of
        only `\sin(a x)` and `\cos(b x)` terms, so these are not implemented.  So,
        for example, you will need to manually convert `\sin^2(x)` into `[1 +
        \cos(2 x)]/2` to properly apply the method of undetermined coefficients on
        it.

        Examples
        ========

        >>> from sympy import log, exp
        >>> from sympy.solvers.ode.ode import _undetermined_coefficients_match
        >>> from sympy.abc import x
        >>> _undetermined_coefficients_match(9*x*exp(x) + exp(-x), x)
        {'test': True, 'trialset': {x*exp(x), exp(-x), exp(x)}}
        >>> _undetermined_coefficients_match(log(x), x)
        {'test': False}

        """
        x = self.ode_problem.sym
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        expr = powsimp(expr, combine='exp')  # exp(x)*exp(2*x + 1) => exp(3*x + 1)
        retdict = {}

        def _test_term(expr, x):
            r"""
            Test if ``expr`` fits the proper form for undetermined coefficients.
            """
            if not expr.has(x):
                return True
            elif expr.is_Add:
                return all(_test_term(i, x) for i in expr.args)
            elif expr.is_Mul:
                if expr.has(sin, cos):
                    foundtrig = False
                    # Make sure that there is only one trig function in the args.
                    # See the docstring.
                    for i in expr.args:
                        if i.has(sin, cos):
                            if foundtrig:
                                return False
                            else:
                                foundtrig = True
                return all(_test_term(i, x) for i in expr.args)
            elif expr.is_Function:
                if expr.func in (sin, cos, exp, sinh, cosh):
                    if expr.args[0].match(a*x + b):
                        return True
                    else:
                        return False
                else:
                    return False
            elif expr.is_Pow and expr.base.is_Symbol and expr.exp.is_Integer and \
                    expr.exp >= 0:
                return True
            elif expr.is_Pow and expr.base.is_number:
                if expr.exp.match(a*x + b):
                    return True
                else:
                    return False
            elif expr.is_Symbol or expr.is_number:
                return True
            else:
                return False

        def _get_trial_set(expr, x, exprs=set()):
            r"""
            Returns a set of trial terms for undetermined coefficients.

            The idea behind undetermined coefficients is that the terms expression
            repeat themselves after a finite number of derivatives, except for the
            coefficients (they are linearly dependent).  So if we collect these,
            we should have the terms of our trial function.
            """
            def _remove_coefficient(expr, x):
                r"""
                Returns the expression without a coefficient.

                Similar to expr.as_independent(x)[1], except it only works
                multiplicatively.
                """
                term = S.One
                if expr.is_Mul:
                    for i in expr.args:
                        if i.has(x):
                            term *= i
                elif expr.has(x):
                    term = expr
                return term

            expr = expand_mul(expr)
            if expr.is_Add:
                for term in expr.args:
                    if _remove_coefficient(term, x) in exprs:
                        pass
                    else:
                        exprs.add(_remove_coefficient(term, x))
                        exprs = exprs.union(_get_trial_set(term, x, exprs))
            else:
                term = _remove_coefficient(expr, x)
                tmpset = exprs.union({term})
                oldset = set()
                while tmpset != oldset:
                    # If you get stuck in this loop, then _test_term is probably
                    # broken
                    oldset = tmpset.copy()
                    expr = expr.diff(x)
                    term = _remove_coefficient(expr, x)
                    if term.is_Add:
                        tmpset = tmpset.union(_get_trial_set(term, x, tmpset))
                    else:
                        tmpset.add(term)
                exprs = tmpset
            return exprs

        def is_homogeneous_solution(term):
            r""" This function checks whether the given trialset contains any root
                of homogenous equation"""
            return expand(sub_func_doit(eq_homogeneous, func, term)).is_zero

        retdict['test'] = _test_term(expr, x)
        if retdict['test']:
            # Try to generate a list of trial solutions that will have the
            # undetermined coefficients. Note that if any of these are not linearly
            # independent with any of the solutions to the homogeneous equation,
            # then they will need to be multiplied by sufficient x to make them so.
            # This function DOES NOT do that (it doesn't even look at the
            # homogeneous equation).
            temp_set = set()
            for i in Add.make_args(expr):
                act = _get_trial_set(i,x)
                if eq_homogeneous is not S.Zero:
                    while any(is_homogeneous_solution(ts) for ts in act):
                        act = {x*ts for ts in act}
                temp_set = temp_set.union(act)

            retdict['trialset'] = temp_set
        return retdict

    def _solve_undetermined_coefficients(self,match):
        r"""
        Helper function for the method of undetermined coefficients.

        See the
        :py:meth:`~sympy.solvers.ode.single.NthConstCoeffUndet`
        docstring for more information on this method.

        The parameter ``match`` should be a dictionary that has the following
        keys:

        ``list``
        A list of solutions to the homogeneous equation.

        ``sol``
        The general solution,.

        ``trialset``
        The set of trial functions as returned by
        ``_undetermined_coefficients_match()['trialset']``.

        """
        eq = self.ode_problem.eq_high_order_free
        x = self.ode_problem.sym
        f = self.ode_problem.func.func
        order = self.ode_problem.order
        r = match
        coeffs = numbered_symbols('a', cls=Dummy)
        coefflist = []
        gensols = r['list']
        gsol = r['sol']
        trialset = r['trialset']
        if len(gensols) != order:
            raise NotImplementedError("Cannot find " + str(order) +
            " solutions to the homogeneous equation necessary to apply" +
            " undetermined coefficients to " + str(eq) +
            " (number of terms != order)")

        trialfunc = 0
        for i in trialset:
            c = next(coeffs)
            coefflist.append(c)
            trialfunc += c*i

        eqs = sub_func_doit(eq, f(x), trialfunc)

        coeffsdict = dict(list(zip(trialset, [0]*(len(trialset) + 1))))

        eqs = _mexpand(eqs)

        for i in Add.make_args(eqs):
            s = separatevars(i, dict=True, symbols=[x])
            if coeffsdict.get(s[x]):
                coeffsdict[s[x]] += s['coeff']
            else:
                coeffsdict[s[x]] = s['coeff']

        coeffvals = solve(list(coeffsdict.values()), coefflist)

        if not coeffvals:
            raise NotImplementedError(
                "Could not solve `%s` using the "
                "method of undetermined coefficients "
                "(unable to solve for coefficients)." % eq)

        psol = trialfunc.subs(coeffvals)

        return Eq(f(x), gsol.rhs + psol)

    def _get_general_solution(self, *, simplify_flag: bool = True):
        eq = self.ode_problem.eq_high_order_free
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        homogen_instance = NthConstCoeffHomogen(SingleODEProblem(eq, f(x), x))
        gensols, collectterms = homogen_instance._get_sols(self.r)
        # A generator of constants
        constants = self.ode_problem.get_numbered_constants(num=len(gensols))
        gsol = Add(*[i*j for (i, j) in zip(constants, gensols)])
        gsol = Eq(f(x), gsol)
        self.r.update({'list': gensols, 'sol': gsol, 'simpliy_flag': simplify_flag})
        gsol = self._solve_undetermined_coefficients(self.r)
        if simplify_flag:
            gsol = homogen_instance._get_simplified_sol([gsol], collectterms)
        return gsol


# Avoid circular import:
from .ode import dsolve, ode_sol_simplicity, odesimp, homogeneous_order
