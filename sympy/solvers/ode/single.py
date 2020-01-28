#
# This is the module for ODE solver classes for single ODEs.
#

from typing import ClassVar, Dict, List, Optional, Type

from sympy.core import Add, S
from sympy.core.exprtools import factor_terms
from sympy.core.expr import Expr
from sympy.core.function import AppliedUndef, Derivative, Function, expand
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Symbol, Dummy, Wild
from sympy.integrals import Integral
from sympy.simplify import simplify, collect
from sympy.utilities import numbered_symbols
from sympy.functions import exp

from sympy.solvers.solvers import solve
from sympy.solvers.deutils import ode_order


class ODEMatchError(NotImplementedError):
    """Raised if a SingleODESolver is asked to solve an ODE it does not match"""
    pass


class SingleODEProblem:
    """Represents an ordinary differential equation (ODE)

    This class is used internally by dsolve and related
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

    def __init__(self, eq, func, sym, eq_orig=None):
        self.eq = eq
        self.eq_orig = eq_orig
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

    def get_dummy_constants(self, num=None):
        return numbered_symbols('C', cls=Dummy, start=1 if num is None else num)


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


class FirstLinear(SingleODESolver):
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

    def _matches(self):
        eq = expand(self.ode_problem.eq)
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        order = ode_order(eq, f(x))
        a = Wild('a', exclude=[f(x)])
        b = Wild('b', exclude=[f(x)])
        c = Wild('c', exclude=[f(x)])
        c1 = Wild('c1', exclude=[x])
        df = f(x).diff(x)

        if order != 1:
            return False

        reduced_eq = None
        if eq.is_Add:
            deriv_coef = eq.coeff(f(x).diff(x, order))
            if deriv_coef not in (1, 0):
                r = deriv_coef.match(a*f(x)**c1)
                if r and r[c1]:
                    den = f(x)**r[c1]
                    reduced_eq = Add(*[arg/den for arg in eq.args])
        if not reduced_eq:
            reduced_eq = eq

        if eq.is_Add:
            ind, dep = reduced_eq.as_independent(f)
        else:
            u = Dummy('u')
            ind, dep = (reduced_eq + u).as_independent(f)
            ind, dep = [tmp.subs(u, 0) for tmp in [ind, dep]]
        coefficients = {a: dep.coeff(df),
                        b: dep.coeff(f(x)),
                        c: ind}
        # double check f[a] since the preconditioning may have failed
        if not coefficients[a].has(f) and not coefficients[b].has(f) and (
                coefficients[a]*df + coefficients[b]*f(x) + coefficients[c]).expand() - reduced_eq == 0:
            coefficients['a'] = a
            coefficients['b'] = b
            coefficients['c'] = c
            self.coeffs = coefficients
            return True

        return False

    def _get_general_solution(self, *, simplify: bool = True):
        x = self.ode_problem.sym
        f = self.ode_problem.func.func
        r = self.coeffs  # a*diff(f(x),x) + b*f(x) + c
        C1 = Dummy('C1')

        t = exp(Integral(r[r['b']]/r[r['a']], x))
        tt = Integral(t*(-r[r['c']]/r[r['a']]), x)
        f = r.get('u', f(x))  # take almost-linear u if present, else f(x)
        return Eq(f, (tt + C1)/t)


class AlmostLinear(FirstLinear):
    r"""
    Solves an almost-linear differential equation.

    The general form of an almost linear differential equation is

    .. math:: f(x) g(y) y + k(x) l(y) + m(x) = 0
                \text{where} l'(y) = g(y)\text{.}

    This can be solved by substituting `l(y) = u(y)`.  Making the given
    substitution reduces it to a linear differential equation of the form `u'
    + P(x) u + Q(x) = 0`.

    The general solution is

        >>> from sympy import Function, dsolve, Eq, pprint
        >>> from sympy.abc import x, y, n
        >>> f, g, k, l = map(Function, ['f', 'g', 'k', 'l'])
        >>> genform = Eq(f(x)*(l(y).diff(y)) + k(x)*l(y) + g(x), 0)
        >>> pprint(genform)
             d
        f(x)*--(l(y)) + g(x) + k(x)*l(y) = 0
             dy
        >>> pprint(dsolve(genform, hint = 'almost_linear'))
               /     //       y*k(x)                \\
               |     ||       ------                ||
               |     ||        f(x)                 ||  -y*k(x)
               |     ||-g(x)*e                      ||  --------
               |     ||--------------  for k(x) != 0||    f(x)
        l(y) = |C1 + |<     k(x)                    ||*e
               |     ||                             ||
               |     ||   -y*g(x)                   ||
               |     ||   --------       otherwise  ||
               |     ||     f(x)                    ||
               \     \\                             //


    See Also
    ========
    :meth:`sympy.solvers.ode.ode.ode_1st_linear`

    Examples
    ========

    >>> from sympy import Function, Derivative, pprint
    >>> from sympy.solvers.ode import dsolve, classify_ode
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> d = f(x).diff(x)
    >>> eq = x*d + x*f(x) + 1
    >>> dsolve(eq, f(x), hint='almost_linear')
    Eq(f(x), (C1 - Ei(x))*exp(-x))
    >>> pprint(dsolve(eq, f(x), hint='almost_linear'))
                         -x
    f(x) = (C1 - Ei(x))*e

    References
    ==========

    - Joel Moses, "Symbolic Integration - The Stormy Decade", Communications
      of the ACM, Volume 14, Number 8, August 1971, pp. 558
    """
    hint = "almost_linear"
    has_integral = True

    def _matches(self):
        ## Almost-linear equation of the form f(x)*g(y)*y' + k(x)*l(y) + m(x) = 0
        eq = expand(self.ode_problem.eq)
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        order = ode_order(eq, f(x))

        if order != 1:
            return False

        df = f(x).diff(x)
        c = Wild('c', exclude=[f(x)])
        d = Wild('d', exclude=[df, f(x).diff(x, 2)])
        e = Wild('e', exclude=[df])

        r = collect(eq, [df, f(x)]).match(e*df + d)
        if r:
            r2 = r.copy()
            r2[c] = S.Zero
            if r2[d].is_Add:
                # Separate the terms having f(x) to r[d] and
                # remaining to r[c]
                no_f, r2[d] = r2[d].as_independent(f(x))
                r2[c] += no_f
            factor = simplify(r2[d].diff(f(x))/r[e])
            if factor and not factor.has(f(x)):
                r2[d] = factor_terms(r2[d])
                u = r2[d].as_independent(f(x), as_Add=False)[1]
                r2.update({'a': e, 'b': d, 'c': c, 'u': u})
                r2[d] /= u
                r2[e] /= u.diff(f(x))
                self.coeffs = r2
                return True

        return False

    def _get_general_solution(self, *, simplify: bool = True):
        return super()._get_general_solution()


class Bernoulli(SingleODESolver):
    r"""
    Solves Bernoulli differential equations.

    These are equations of the form

    .. math:: dy/dx + P(x) y = Q(x) y^n\text{, }n \ne 1`\text{.}

    The substitution `w = 1/y^{1-n}` will transform an equation of this form
    into one that is linear (see the docstring of
    :py:meth:`~sympy.solvers.ode.ode.ode_1st_linear`).  The general solution is::

        >>> from sympy import Function, dsolve, Eq, pprint
        >>> from sympy.abc import x, n
        >>> f, P, Q = map(Function, ['f', 'P', 'Q'])
        >>> genform = Eq(f(x).diff(x) + P(x)*f(x), Q(x)*f(x)**n)
        >>> pprint(genform)
                    d                n
        P(x)*f(x) + --(f(x)) = Q(x)*f (x)
                    dx
        >>> pprint(dsolve(genform, f(x), hint='Bernoulli_Integral'), num_columns=100)
                                                                                      1
                                                                                    -----
                                                                                    1 - n
               //               /                            \                     \
               ||              |                             |                     |
               ||              |                  /          |             /       |
               ||              |                 |           |            |        |
               ||              |        (1 - n)* | P(x) dx   |  -(1 - n)* | P(x) dx|
               ||              |                 |           |            |        |
               ||              |                /            |           /         |
        f(x) = ||C1 + (n - 1)* | -Q(x)*e                   dx|*e                   |
               ||              |                             |                     |
               \\              /                            /                     /


    Note that the equation is separable when `n = 1` (see the docstring of
    :py:meth:`~sympy.solvers.ode.ode.ode_separable`).

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
    f(x) = -------------------
             /     log(x)   1\
           x*|C1 + ------ + -|
             \       x      x/

    References
    ==========

    - https://en.wikipedia.org/wiki/Bernoulli_differential_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 95

    # indirect doctest

    """
    hint = "Bernoulli"
    has_integral = True

    def _matches(self):
        eq = expand(self.ode_problem.eq)
        f = self.ode_problem.func.func
        x = self.ode_problem.sym
        order = ode_order(eq, f(x))
        df = f(x).diff(x)

        if order != 1:
            return False

        a = Wild('a', exclude=[f(x)])
        b = Wild('b', exclude=[f(x)])
        c = Wild('c', exclude=[f(x)])
        n = Wild('n', exclude=[x, f(x), df])
        c1 = Wild('c1', exclude=[x])

        reduced_eq = None
        if eq.is_Add:
            deriv_coef = eq.coeff(f(x).diff(x, order))
            if deriv_coef not in (1, 0):
                r = deriv_coef.match(a*f(x)**c1)
                if r and r[c1]:
                    den = f(x)**r[c1]
                    reduced_eq = Add(*[arg/den for arg in eq.args])
        if not reduced_eq:
            reduced_eq = eq

        r = collect(
            reduced_eq, f(x), exact=True).match(a*df + b*f(x) + c*f(x)**n)
        if r and r[c] != 0 and r[n] != 1:  # See issue 4676
            r['a'] = a
            r['b'] = b
            r['c'] = c
            r['n'] = n
            self.coeffs = r
            return True
        return False

    def _get_general_solution(self, *, simplify: bool = True):
        x = self.ode_problem.sym
        f = self.ode_problem.func.func
        r = self.coeffs  # a*diff(f(x),x) + b*f(x) + c*f(x)**n, n != 1
        C1 = Dummy('C1')
        t = exp((1 - r[r['n']])*Integral(r[r['b']]/r[r['a']], x))
        tt = (r[r['n']] - 1)*Integral(t*r[r['c']]/r[r['a']], x)
        return Eq(f(x), ((tt + C1)/t)**(1/(1 - r[r['n']])))
