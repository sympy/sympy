"""Solvers of systems of polynomial equations. """

from sympy.polys import Poly, groebner, roots
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.polys.polyerrors import ComputationFailed
from sympy.utilities import any, all, postfixes
from sympy.utilities.iterables import minkey
from sympy.simplify import rcollect
from sympy.core import S

class SolveFailed(Exception):
    """Raised when solver's conditions weren't met. """

def solve_poly_system(seq, *gens, **args):
    """
    Solve a system of polynomial equations.

    Example
    =======

    >>> from sympy import solve_poly_system
    >>> from sympy.abc import x, y

    >>> solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y)
    [(0, 0), (2, -2**(1/2)), (2, 2**(1/2))]

    """
    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('solve_poly_system', len(seq), exc)

    if len(polys) == len(opt.gens) == 2:
        f, g = polys

        a, b = f.degree_list()
        c, d = g.degree_list()

        if a <= 2 and b <= 2 and c <= 2 and d <= 2:
            try:
                return solve_biquadratic(f, g, opt)
            except SolveFailed:
                pass

    return solve_generic(polys, opt)

def solve_biquadratic(f, g, opt):
    """Solve a system of two bivariate quadratic polynomial equations. """
    G = groebner([f, g])

    if len(G) == 1 and G[0].is_ground:
        return None

    if len(G) != 2:
        raise SolveFailed

    p, q = G
    x, y = opt.gens

    p = Poly(p, x, expand=False)
    q = q.ltrim(-1)

    p_roots = [ rcollect(expr, y) for expr in roots(p).keys() ]
    q_roots = roots(q).keys()

    solutions = []

    for q_root in q_roots:
        for p_root in p_roots:
            solution = (p_root.subs(y, q_root), q_root)
            solutions.append(solution)

    return sorted(solutions)

def solve_generic(polys, opt):
    """
    Solve a generic system of polynomial equations.

    Returns all possible solutions over C[x_1, x_2, ..., x_m] of a
    set F = { f_1, f_2, ..., f_n } of polynomial equations,  using
    Groebner basis approach. For now only zero-dimensional systems
    are supported, which means F can have at most a finite number
    of solutions.

    The algorithm works by the fact that, supposing G is the basis
    of F with respect to an elimination order  (here lexicographic
    order is used), G and F generate the same ideal, they have the
    same set of solutions. By the elimination property,  if G is a
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

    References
    ==========

    .. [Buchberger01] B. Buchberger, Groebner Bases: A Short
    Introduction for Systems Theorists, In: R. Moreno-Diaz,
    B. Buchberger, J.L. Freire, Proceedings of EUROCAST'01,
    February, 2001

    .. [Cox97] D. Cox, J. Little, D. O'Shea, Ideals, Varieties
    and Algorithms, Springer, Second Edition, 1997, pp. 112

    """
    def is_univariate(f):
        """Returns True if 'f' is univariate in its last variable. """
        for monom in f.monoms():
            if any(m > 0 for m in monom[:-1]):
                return False

        return True

    def subs_root(f, gen, zero):
        """Replace generator with a root so that the result is nice. """
        p = f.as_expr({gen: zero})

        if f.degree(gen) >= 2:
            p = p.expand(deep=False)

        return p

    def solve_reduced_system(system, gens, entry=False):
        """Recursively solves reduced polynomial systems. """
        if len(system) == len(gens) == 1:
            zeros = roots(system[0], gens[-1]).keys()
            return [ (zero,) for zero in zeros ]

        basis = groebner(system, gens, polys=True)

        if len(basis) == 1 and basis[0].is_ground:
            if not entry:
                return []
            else:
                return None

        univariate = filter(is_univariate, basis)

        if len(univariate) == 1:
            f = univariate.pop()
        else:
            raise NotImplementedError("only zero-dimensional systems supported (finite number of solutions)")

        gens = f.gens
        gen = gens[-1]

        zeros = roots(f.ltrim(gen)).keys()

        if not zeros:
            return []

        if len(basis) == 1:
            return [ (zero,) for zero in zeros ]

        solutions = []

        for zero in zeros:
            new_system = []
            new_gens = gens[:-1]

            for b in basis[:-1]:
                eq = subs_root(b, gen, zero)

                if eq is not S.Zero:
                    new_system.append(eq)

            for solution in solve_reduced_system(new_system, new_gens):
                solutions.append(solution + (zero,))

        return solutions

    result = solve_reduced_system(polys, opt.gens, entry=True)

    if result is not None:
        return sorted(result)
    else:
        return None

def solve_triangulated(polys, *gens, **args):
    """
    Solve a polynomial system using Gianni-Kalkbrenner algorithm.

    The algorithm proceeds by computing one Groebner basis in the ground
    domain and then by iteratively computing polynomial factorizations in
    appropriately constructed algebraic extensions of the ground domain.

    Example
    =======

    >>> from sympy.solvers.polysys import solve_triangulated
    >>> from sympy.abc import x, y, z

    >>> F = [x**2 + y + z - 1, x + y**2 + z - 1, x + y + z**2 - 1]

    >>> solve_triangulated(F, x, y, z)
    [(0, 0, 1), (0, 1, 0), (1, 0, 0)]

    References
    ==========

    .. [Gianni89] Patrizia Gianni, Teo Mora, Algebraic Solution of System of
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
    solutions = set([])

    for zero in zeros:
        solutions.add(((zero,), dom))

    var_seq = reversed(gens[:-1])
    vars_seq = postfixes(gens[1:])

    for var, vars in zip(var_seq, vars_seq):
        _solutions = set([])

        for values, dom in solutions:
            H, mapping = [], zip(vars, values)

            for g in G:
                _vars = (var,) + vars

                if g.has_only_gens(*_vars) and g.degree(var) != 0:
                    h = g.ltrim(var).eval(mapping)

                    if g.degree(var) == h.degree():
                        H.append(h)

            p = minkey(H, key=lambda h: h.degree())
            zeros = p.ground_roots()

            for zero in zeros:
                if not zero.is_Rational:
                    dom_zero = dom.algebraic_field(zero)
                else:
                    dom_zero = dom

                _solutions.add(((zero,) + values, dom_zero))

        solutions = _solutions

    solutions = list(solutions)

    for i, (solution, _) in enumerate(solutions):
        solutions[i] = solution

    return sorted(solutions)
