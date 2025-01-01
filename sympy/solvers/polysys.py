"""Solvers of systems of polynomial equations. """
import itertools
from collections import defaultdict

from sympy import Dummy
from sympy.core import S
from sympy.core.sorting import default_sort_key
from sympy.polys import Poly, groebner, roots
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.polys.polyerrors import (
    ComputationFailed,
    PolificationFailed,
    CoercionFailed,
    GeneratorsNeeded,
    DomainError
)
from sympy.simplify import rcollect
from sympy.utilities import postfixes
from sympy.utilities.misc import filldedent
from sympy.utilities.iterables import cartes
from sympy.logic.boolalg import Or, And
from sympy.core.relational import Eq


class SolveFailed(Exception):
    """Raised when solver's conditions were not met. """


def solve_poly_system(seq, *gens, strict=False, **args):
    """
    Return a list of solutions for the system of polynomial equations
    or else None.

    Parameters
    ==========

    seq: a list/tuple/set
        Listing all the equations that are needed to be solved
    gens: generators
        generators of the equations in seq for which we want the
        solutions
    strict: a boolean (default is False)
        if strict is True, NotImplementedError will be raised if
        the solution is known to be incomplete (which can occur if
        not all solutions are expressible in radicals)
    args: Keyword arguments
        Special options for solving the equations.


    Returns
    =======

    List[Tuple]
        a list of tuples with elements being solutions for the
        symbols in the order they were passed as gens
    None
        None is returned when the computed basis contains only the ground.

    Examples
    ========

    >>> from sympy import solve_poly_system
    >>> from sympy.abc import x, y

    >>> solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y)
    [(0, 0), (2, -sqrt(2)), (2, sqrt(2))]

    >>> solve_poly_system([x**5 - x + y**3, y**2 - 1], x, y, strict=True)
    Traceback (most recent call last):
    ...
    UnsolvableFactorError

    """
    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('solve_poly_system', len(seq), exc)

    if len(polys) == len(opt.gens) == 2:
        f, g = polys

        if all(i <= 2 for i in f.degree_list() + g.degree_list()):
            try:
                return solve_biquadratic(f, g, opt)
            except SolveFailed:
                pass

    return solve_generic(polys, opt, strict=strict)


def solve_biquadratic(f, g, opt):
    """Solve a system of two bivariate quadratic polynomial equations.

    Parameters
    ==========

    f: a single Expr or Poly
        First equation
    g: a single Expr or Poly
        Second Equation
    opt: an Options object
        For specifying keyword arguments and generators

    Returns
    =======

    List[Tuple]
        a list of tuples with elements being solutions for the
        symbols in the order they were passed as gens
    None
        None is returned when the computed basis contains only the ground.

    Examples
    ========

    >>> from sympy import Options, Poly
    >>> from sympy.abc import x, y
    >>> from sympy.solvers.polysys import solve_biquadratic
    >>> NewOption = Options((x, y), {'domain': 'ZZ'})

    >>> a = Poly(y**2 - 4 + x, y, x, domain='ZZ')
    >>> b = Poly(y*2 + 3*x - 7, y, x, domain='ZZ')
    >>> solve_biquadratic(a, b, NewOption)
    [(1/3, 3), (41/27, 11/9)]

    >>> a = Poly(y + x**2 - 3, y, x, domain='ZZ')
    >>> b = Poly(-y + x - 4, y, x, domain='ZZ')
    >>> solve_biquadratic(a, b, NewOption)
    [(7/2 - sqrt(29)/2, -sqrt(29)/2 - 1/2), (sqrt(29)/2 + 7/2, -1/2 + \
      sqrt(29)/2)]
    """
    G = groebner([f, g])

    if len(G) == 1 and G[0].is_ground:
        return None

    if len(G) != 2:
        raise SolveFailed

    x, y = opt.gens
    p, q = G
    if not p.gcd(q).is_ground:
        # not 0-dimensional
        raise SolveFailed

    p = Poly(p, x, expand=False)
    p_roots = [rcollect(expr, y) for expr in roots(p).keys()]

    q = q.ltrim(-1)
    q_roots = list(roots(q).keys())

    solutions = [(p_root.subs(y, q_root), q_root) for q_root, p_root in
                 itertools.product(q_roots, p_roots)]

    return sorted(solutions, key=default_sort_key)


def solve_generic(polys, opt, strict=False):
    """
    Solve a generic system of polynomial equations.

    Returns all possible solutions over C[x_1, x_2, ..., x_m] of a
    set F = { f_1, f_2, ..., f_n } of polynomial equations, using
    Groebner basis approach. For now only zero-dimensional systems
    are supported, which means F can have at most a finite number
    of solutions. If the basis contains only the ground, None is
    returned.

    The algorithm works by the fact that, supposing G is the basis
    of F with respect to an elimination order (here lexicographic
    order is used), G and F generate the same ideal, they have the
    same set of solutions. By the elimination property, if G is a
    reduced, zero-dimensional Groebner basis, then there exists an
    univariate polynomial in G (in its last variable). This can be
    solved by computing its roots. Substituting all computed roots
    for the last (eliminated) variable in other elements of G, new
    polynomial system is generated. Applying the above procedure
    recursively, a finite number of solutions can be found.

    The ability of finding all solutions by this procedure depends
    on the root finding algorithms. If no solutions were found, it
    means only that roots() failed, but the system is solvable. To
    overcome this difficulty use numerical algorithms instead.

    Parameters
    ==========

    polys: a list/tuple/set
        Listing all the polynomial equations that are needed to be solved
    opt: an Options object
        For specifying keyword arguments and generators
    strict: a boolean
        If strict is True, NotImplementedError will be raised if the solution
        is known to be incomplete

    Returns
    =======

    List[Tuple]
        a list of tuples with elements being solutions for the
        symbols in the order they were passed as gens
    None
        None is returned when the computed basis contains only the ground.

    References
    ==========

    .. [Buchberger01] B. Buchberger, Groebner Bases: A Short
    Introduction for Systems Theorists, In: R. Moreno-Diaz,
    B. Buchberger, J.L. Freire, Proceedings of EUROCAST'01,
    February, 2001

    .. [Cox97] D. Cox, J. Little, D. O'Shea, Ideals, Varieties
    and Algorithms, Springer, Second Edition, 1997, pp. 112

    Raises
    ========

    NotImplementedError
        If the system is not zero-dimensional (does not have a finite
        number of solutions)

    UnsolvableFactorError
        If ``strict`` is True and not all solution components are
        expressible in radicals

    Examples
    ========

    >>> from sympy import Poly, Options
    >>> from sympy.solvers.polysys import solve_generic
    >>> from sympy.abc import x, y
    >>> NewOption = Options((x, y), {'domain': 'ZZ'})

    >>> a = Poly(x - y + 5, x, y, domain='ZZ')
    >>> b = Poly(x + y - 3, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [(-1, 4)]

    >>> a = Poly(x - 2*y + 5, x, y, domain='ZZ')
    >>> b = Poly(2*x - y - 3, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [(11/3, 13/3)]

    >>> a = Poly(x**2 + y, x, y, domain='ZZ')
    >>> b = Poly(x + y*4, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [(0, 0), (1/4, -1/16)]

    >>> a = Poly(x**5 - x + y**3, x, y, domain='ZZ')
    >>> b = Poly(y**2 - 1, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption, strict=True)
    Traceback (most recent call last):
    ...
    UnsolvableFactorError

    """
    def _is_univariate(f):
        """Returns True if 'f' is univariate in its last variable. """
        for monom in f.monoms():
            if any(monom[:-1]):
                return False

        return True

    def _subs_root(f, gen, zero):
        """Replace generator with a root so that the result is nice. """
        p = f.as_expr({gen: zero})

        if f.degree(gen) >= 2:
            p = p.expand(deep=False)

        return p

    def _solve_reduced_system(system, gens, entry=False):
        """Recursively solves reduced polynomial systems. """
        if len(system) == len(gens) == 1:
            # the below line will produce UnsolvableFactorError if
            # strict=True and the solution from `roots` is incomplete
            zeros = list(roots(system[0], gens[-1], strict=strict).keys())
            return [(zero,) for zero in zeros]

        basis = groebner(system, gens, polys=True)

        if len(basis) == 1 and basis[0].is_ground:
            if not entry:
                return []
            else:
                return None

        univariate = list(filter(_is_univariate, basis))

        if len(basis) < len(gens):
            raise NotImplementedError(filldedent('''
                only zero-dimensional systems supported
                (finite number of solutions)
                '''))

        if len(univariate) == 1:
            f = univariate.pop()
        else:
            raise NotImplementedError(filldedent('''
                only zero-dimensional systems supported
                (finite number of solutions)
                '''))

        gens = f.gens
        gen = gens[-1]

        # the below line will produce UnsolvableFactorError if
        # strict=True and the solution from `roots` is incomplete
        zeros = list(roots(f.ltrim(gen), strict=strict).keys())

        if not zeros:
            return []

        if len(basis) == 1:
            return [(zero,) for zero in zeros]

        solutions = []

        for zero in zeros:
            new_system = []
            new_gens = gens[:-1]

            for b in basis[:-1]:
                eq = _subs_root(b, gen, zero)

                if eq is not S.Zero:
                    new_system.append(eq)

            for solution in _solve_reduced_system(new_system, new_gens):
                solutions.append(solution + (zero,))

        if solutions and len(solutions[0]) != len(gens):
            raise NotImplementedError(filldedent('''
                only zero-dimensional systems supported
                (finite number of solutions)
                '''))
        return solutions

    try:
        result = _solve_reduced_system(polys, opt.gens, entry=True)
    except CoercionFailed:
        raise NotImplementedError

    if result is not None:
        return sorted(result, key=default_sort_key)


def solve_triangulated(polys, *gens, **args):
    """
    Solve a polynomial system using Gianni-Kalkbrenner algorithm.

    The algorithm proceeds by computing one Groebner basis in the ground
    domain and then by iteratively computing polynomial factorizations in
    appropriately constructed algebraic extensions of the ground domain.

    Parameters
    ==========

    polys: a list/tuple/set
        Listing all the equations that are needed to be solved
    gens: generators
        generators of the equations in polys for which we want the
        solutions
    args: Keyword arguments
        Special options for solving the equations

    Returns
    =======

    List[Tuple]
        A List of tuples. Solutions for symbols that satisfy the
        equations listed in polys

    Examples
    ========

    >>> from sympy import solve_triangulated
    >>> from sympy.abc import x, y, z

    >>> F = [x**2 + y + z - 1, x + y**2 + z - 1, x + y + z**2 - 1]

    >>> solve_triangulated(F, x, y, z)
    [(0, 0, 1), (0, 1, 0), (1, 0, 0)]

    References
    ==========

    1. Patrizia Gianni, Teo Mora, Algebraic Solution of System of
    Polynomial Equations using Groebner Bases, AAECC-5 on Applied Algebra,
    Algebraic Algorithms and Error-Correcting Codes, LNCS 356 247--257, 1989

    """
    G = groebner(polys, gens, polys=True)
    G = list(reversed(G))

    domain = args.get('domain')

    if domain is not None:
        for i, g in enumerate(G):
            G[i] = g.set_domain(domain)

    f, G = G[0].ltrim(-1), G[1:]
    dom = f.get_domain()

    zeros = f.ground_roots()
    solutions = {((zero,), dom) for zero in zeros}

    var_seq = reversed(gens[:-1])
    vars_seq = postfixes(gens[1:])

    for var, vars in zip(var_seq, vars_seq):
        _solutions = set()

        for values, dom in solutions:
            H, mapping = [], list(zip(vars, values))

            for g in G:
                _vars = (var,) + vars

                if g.has_only_gens(*_vars) and g.degree(var) != 0:
                    h = g.ltrim(var).eval(dict(mapping))

                    if g.degree(var) == h.degree():
                        H.append(h)

            p = min(H, key=lambda h: h.degree())
            zeros = p.ground_roots()

            for zero in zeros:
                if not zero.is_Rational:
                    dom_zero = dom.algebraic_field(zero)
                else:
                    dom_zero = dom

                _solutions.add(((zero,) + values, dom_zero))

        solutions = _solutions
    return sorted((s for s, _ in solutions), key=default_sort_key)


def factor_system(eqs, *gens, **kwargs):
    """
    Factorizes a system of polynomial equations into
    irreducible subsystems.

    Parameters
    ==========

    eqs : list
        List of expressions to be factored.
        Each expression is assumed to be equal to zero.

    *gens : Symbol or sequence of Symbols, optional
        Generator(s) of the polynomial ring.
        If not provided, all free symbols will be used.

    **kwargs : dict, optional
        Same optional arguments taken by ``factor``

    Returns
    =======

    tuple
        A pair (systems, condition) where:

        * systems is a list of lists of expressions,
          where each sublist when solved gives one
          component of the solution

        * condition is a Boolean expression that must
          be satisfied for the system to be solvable

    Examples
    ========

    >>> from sympy.solvers.polysys import factor_system
    >>> from sympy.abc import x, y, a, b, c

    A simple system with multiple solutions:

    >>> systems, cond = factor_system([x**2 - 1, y - 1])
    >>> systems
    [[x - 1, y - 1], [x + 1, y - 1]]
    >>> cond
    True

    A system with no solution:

    >>> systems, cond = factor_system([x, 1])
    >>> systems
    []
    >>> cond
    True

    >>> systems, cond = factor_system([a*x*(x-1), b*y, c], [x, y])
    >>> systems
    [[x, y], [x - 1, y]]
    >>> cond
    Eq(c, 0)

    The above cond signifies that the symbol c must be equal to
    zero for the system to be solvable.

    A system with infinite solutions:

    >>> systems, cond = factor_system([x - x, (x + 1)**2 - (x**2 + 2*x + 1)])
    >>> systems
    [[]]
    >>> cond
    True

    The solution set of the original system represented
    by eqs is the union of the solution sets of the
    factorized systems.

    A return of a list containing an empty list [[]]
    in the tuple means any value of the symbol(s)
    is a solution whereas a return of an empty list []
    in the tuple means no solutions for the system exists.

    See Also
    ========

    factor_system_cond : Returns both factors and conditions
    factor_system_bool : Returns a Boolean combination representing the solution
    sympy.polys.polytools.factor : Factors a polynomial into irreducible factors
                                   over the rational numbers
    """

    systems, conds = factor_system_cond(eqs, *gens, **kwargs)
    systems = [[f for f, c in system] for system in systems]
    return systems, And(*[Eq(c, 0) for c in conds]) if conds else True


def factor_system_bool(eqs, *gens, **kwargs):
    """
    Factorizes a system of polynomial equations into irreducible DNF.

    The system of expressions(eqs) is taken and a Boolean combination
    of equations is returned that represents the same solution set.
    The result is in disjunctive normal form (OR of ANDs).

    Parameters
    ==========

    eqs : list
       List of expressions to be factored.
       Each expression is assumed to be equal to zero.

    *gens : Symbol or sequence of Symbols, optional
       Generator(s) of the polynomial ring.
       If not provided, all free symbols will be used.

    **kwargs : dict, optional
       Optional keyword arguments


    Returns
    =======

    Boolean:
       A Boolean combination of equations. The result is typically in
       the form of a conjunction (AND) of a disjunctive normal form
       with additional conditions.

    Examples
    ========

    >>> from sympy.solvers.polysys import factor_system_bool
    >>> from sympy.abc import x, y, a, b, c
    >>> factor_system_bool([x**2 - 1])
    Eq(x - 1, 0) | Eq(x + 1, 0)

    >>> factor_system_bool([x**2 - 1, y - 1])
    (Eq(x - 1, 0) & Eq(y - 1, 0)) | (Eq(x + 1, 0) & Eq(y - 1, 0))

    >>> eqs = [a * (x - 1), b]
    >>> factor_system_bool(eqs, [x])
    Eq(b, 0) & (Eq(a, 0) | Eq(x - 1, 0))

    >>> eqs = [a * x ** 2 - a, b * (x + 1), c]
    >>> factor_system_bool(eqs, [x])
    Eq(c, 0) & (Eq(x + 1, 0) | (Eq(a, 0) & Eq(b, 0)))

    >>> factor_system_bool([x**2 + 2*x + 1 - (x + 1)**2])
    True

    The result is logically equivalent to the system of equations
    i.e. eqs. The function returns ``True`` when all values of
    the symbol(s) is a solution and ``False`` when the system
    cannot be solved.

    See Also
    ========

    factor_system : Returns factors and solvability condition separately
    factor_system_cond : Returns both factors and conditions

    """

    systems, conds = factor_system_cond(eqs, *gens, **kwargs)
    sys = Or(*[And(*[_eq2bool(eq) for eq in sys]) for sys in systems])
    if conds:
        sys &= And(*[Eq(c, 0) for c in conds])
    return sys


def _eq2bool(eq):
    """Convert a (factor, conditions) pair to a Boolean equation."""
    f, cs = eq
    b = Eq(f, 0)
    if cs:
        b |= And(*[Eq(c, 0) for c in cs])
    return b


def factor_system_cond(eqs, *gens, **kwargs):
    """
    Factorizes a polynomial system into irreducible components and returns
    factors with their associated conditions.

    Parameters
    ==========

    eqs : list
        List of expressions to be factored.
        Each expression is assumed to be equal to zero.

    *gens : Symbol or sequence of Symbols, optional
        Generator(s) of the polynomial ring.
        If not provided, all free symbols will be used.

    **kwargs : dict, optional
        Optional keyword arguments.

    Returns
    =======

    tuple
        A pair (systems, conditions) where:

        * systems is a list of lists of (factor, conditions) pairs,
          where factor is an expression and conditions is a list
          of expressions that must be zero for the factor to be relevant

        * conditions is a list of expressions that must be zero
          for any solution to exist

    Examples
    ========

    >>> from sympy.solvers.polysys import factor_system_cond
    >>> from sympy.abc import x, y, a, b, c

    >>> systems, conds = factor_system_cond([x**2 - 4, a*y, b], [x, y])
    >>> systems
    [[(x - 2, []), (y, [a])], [(x + 2, []), (y, [a])]]
    >>> conds
    [b]

    >>> systems, conds = factor_system_cond([a*x*(x-1), b*y, c], [x, y])
    >>> systems
    [[(x, [a]), (y, [b])], [(x - 1, [a]), (y, [b])]]
    >>> conds
    [c]

    In the return of the (systems, conds) pair, systems is [[]]
    when the system implies tautology. Eg. x - x = 0,
    and [] when the system of equations is unatisfiable. Eg. 1 = 0

    See Also
    ========

    factor_system : Returns factors and solvability condition separately
    factor_system_bool : Returns a Boolean combination representing the solution
    sympy.polys.polytools.factor : Factors a polynomial into irreducible factors
                                   over the rational numbers
    """
    try:
        polys, opts = parallel_poly_from_expr(eqs, *gens, **kwargs)
        only_numbers = False
    except (GeneratorsNeeded, PolificationFailed):
        _u = Dummy('u')
        polys, opts = parallel_poly_from_expr(eqs, [_u], **kwargs)
        assert opts['domain'].is_Numerical
        only_numbers = True

    if only_numbers:
        if all(p == 0 for p in polys):
            systems = [[]]
            conditions = []
        else:
            systems = []
            conditions = []
    else:
        systems, conditions = factor_system_poly(polys)
        systems = [[(f.as_expr(), [c.as_expr() for c in cs]) for f, cs in system] for system in systems]
        conditions = [c.as_expr() for c in conditions]

    return systems, conditions


def factor_system_poly(polys):
    """
    Factors a system of polynomials into irreducible factors with conditions.
    Core implementation that works directly with Poly instances.

    Parameters
    ==========

    polys : list
        A list of Poly instances to be factored.

    Returns
    =======

    tuple
        A pair (systems, constant_conds) where:

        * systems is a list of lists of (factor, conditions) pairs, with each
          factor being a Poly instance and conditions being a tuple of Poly
          instances

        * constant_conds is a list of Poly instances representing conditions that
          must be satisfied for any solution to exist

    Examples
    ========

    >>> from sympy import symbols, Poly, ZZ
    >>> from sympy.solvers.polysys import factor_system_poly
    >>> a, b, c, x = symbols('a b c x')
    >>> p1 = Poly((a - 1)*(x - 2), x, domain=ZZ[a,b,c])
    >>> p2 = Poly((b - 3)*(x - 2), x, domain=ZZ[a,b,c])
    >>> p3 = Poly(c, x, domain=ZZ[a,b,c])

    The equation to be solved for x is ``x - 2 = 0`` provided either
    of the two conditions on the parameters ``a`` and ``b`` is nonzero
    and the constant parameter ``c`` should be zero.

    >>> systems, constant_conds = factor_system_poly([p1, p2, p3])
    >>> systems
    [[(Poly(x - 2, x, domain='ZZ[a,b,c]'), (Poly(a - 1, x, domain='ZZ[a,b,c]'), Poly(b - 3, x, domain='ZZ[a,b,c]')))]]

    >>> constant_conds
    [Poly(c, x, domain='ZZ[a,b,c]')]

    This is the core routine used by higher-level functions factor_system_cond,
    factor_system, and factor_system_bool. It Returns empty systems list
    ([]) when no solution exists and [[]] (list containing empty list) when any
    value is a solution

    See Also
    ========

    factor_system : Returns factors and solvability condition separately
    factor_system_bool : Returns a Boolean combination representing the solution
    factor_system_cond : Returns both factors and conditions
    sympy.polys.polytools.factor : Factors a polynomial into irreducible factors
                                   over the rational numbers
    """
    if not all(isinstance(poly, Poly) for poly in polys):
        raise TypeError("polys should be a list of Poly instances")
    if not polys:
        return [[]], []
    else:
        base_domain = polys[0].domain
        base_gens = polys[0].gens
        if not all(poly.domain == base_domain and poly.gens == base_gens for poly in polys[1:]):
            raise DomainError("All polynomials must have the same domain and generators")

    constant_eqs = []
    eqs_factors = []
    conds_factor = defaultdict(list)

    for poly in polys:
        constant, factors_mult = poly.factor_list()

        factors = [f for f, m in factors_mult]

        if factors:
            eqs_factors.append(factors)
            if constant.is_zero is not False:
                constp = constant.as_poly(base_gens, domain=base_domain)
                for f in factors:
                    if constp not in conds_factor[f]:
                        conds_factor[f].append(constp)
        elif constant.is_zero is True:
            pass
        elif constant.is_zero is False:
            return ([], [])
        else:
            constp = constant.as_poly(base_gens, domain=base_domain)
            constant_eqs.append(constp)

    fac2conds = {f: tuple(conds) for f, conds in conds_factor.items()}

    cnf = [[(f, fac2conds.get(f, ())) for f in eq] for eq in eqs_factors]

    dnf = _cnf2dnf(cnf)

    return dnf, constant_eqs


def _has_common_variables(poly1, poly2):
    """Helper function to check if two polynomials have nonzero degree in any common variables"""
    return any(d1 and d2 for d1, d2 in zip(poly1.degree_list(), poly2.degree_list()))


def _cnf2dnf(eqs):
    """
    Given a list of lists of (factor, conditions) pairs from the factorization of a
    polynomial system, returns the minimal DNF sufficient to
    satisfy the CNF. Only includes terms necessary for
    satisfying the CNF, omitting redundant factors, hence different from a simple
    mechanical rewrite of CNF to DNF.

    The input is a list of lists of (Poly, tuple) pairs, where each Poly represents a
    factor and the tuple contains conditions on that factor..
    """
    # remove equations that are independent of all others
    # and generate the Cartesian product for them
    eqs = [list(i) for i in eqs]
    indep = set()

    for i in range(len(eqs)):
        is_independent = True
        for j in range(len(eqs)):
            if i != j:
                for expr1 in eqs[i]:
                    for expr2 in eqs[j]:
                        if _has_common_variables(expr1[0], expr2[0]):
                            is_independent = False
                            break
                    if not is_independent:
                        break
            if not is_independent:
                break
        if is_independent:
            indep.add(i)

    if not indep:
        return list(_dnf(eqs))

    indep_eqs = []
    dep_eqs = []
    for i, e in enumerate(eqs):
        if i in indep:
            indep_eqs.append(e)
        else:
            dep_eqs.append(e)

    result = []
    for t in cartes(*indep_eqs):
        t_list = list(t)
        for d in _dnf(dep_eqs):
            result.append(t_list + d)
    return result


def _dnf(eqs):
    # helper for _cnf2dnf that recursively enumerates the
    # minimal dnf args that satisfy the cnf expression;
    if not eqs:
        return [[]]
    elif len(eqs) == 1:
        return [[f] for f in eqs[0]]
    else:
        f = eqs[0][0]
        eqs_f_zero = [eq for eq in eqs if f not in eq]
        dnf_f = _dnf(eqs_f_zero)

        f_free_eqs = [[x for x in eq if x != f] for eq in eqs]
        if not all(f_free_eqs):
            dnf = [[f] + s for s in dnf_f]
        else:
            eqs_f_nonzero = list(filter(None, f_free_eqs))
            dnf_no_f = _dnf(eqs_f_nonzero)
            dnf = dnf_no_f + [[f] + s for s in dnf_f if s not in dnf_no_f]

        return dnf
