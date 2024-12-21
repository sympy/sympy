"""Solvers of systems of polynomial equations. """
import itertools
from collections import defaultdict

from sympy import Dummy
from sympy.core import S
from sympy.core.sorting import default_sort_key
from sympy.polys import Poly, groebner, roots
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.polys.polyerrors import (ComputationFailed,
                                    PolificationFailed, CoercionFailed, GeneratorsNeeded, DomainError)
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


def factor_poly_system(eqs, *gens, **kwargs):

    systems, conds = factor_system_cond(eqs, *gens, **kwargs)
    systems = [[f.as_expr() for f, c in system] for system in systems]
    return systems, And(*conds)


def factor_system_bool(eqs, *gens, **kwargs):

    systems, conds = factor_system_cond(eqs, *gens, **kwargs)
    sys = Or(*[And(*[_eq2bool(eq) for eq in sys]) for sys in systems])
    if conds:
        sys &= And(*[Eq(c.as_expr(), 0) for c in conds])
    return sys


def _eq2bool(eq):
    """Convert a (factor, conditions) pair to a Boolean equation."""
    f, cs = eq
    b = Eq(f.as_expr(), 0)
    if cs:
        b |= And(*[Eq(c.as_expr(), 0) for c in cs])
    return b


def factor_system_cond(eqs, *gens, **kwargs):

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

    return systems, conditions


def factor_system_poly(polys):

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


def _cnf2dnf(eqs):
    # remove equations that are independent of all others
    # and generate the Cartesian product for them
    eqs = [list(i) for i in eqs]
    indep = []
    for i in range(len(eqs)):
        for j in range(len(eqs)):
            if i != j and any(x[0] == y[0] for x in eqs[i] for y in eqs[j]):
                break
        else:
            indep.append(i)

    if not indep:
        return list(_dnf(eqs))
    else:
        result = []
        indep_eqs = [eqs.pop(i) for i in indep[::-1]]
        for t in cartes(*indep_eqs):
            t_list = list(t)
            for d in _dnf(eqs):
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
