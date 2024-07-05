"""Solvers of systems of polynomial equations. """
import itertools

from sympy.core import S, N
from sympy.core.sorting import default_sort_key
from sympy.core.function import diff
from sympy.polys import Poly, groebner, roots, real_roots, nroots
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.polys.polytools import LT, LC, degree, subresultants
from sympy.polys.rootoftools import ComplexRootOf
from sympy.polys.polyerrors import (ComputationFailed,
    PolificationFailed, CoercionFailed)
from sympy.simplify import rcollect
from sympy.utilities import postfixes
from sympy.utilities.misc import filldedent


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


def solve_poly_system_cad(seq, gens, return_one_sample=True):
    """
    Solves a system of polynomial inequalities/equalities via
    cylindrical algebraic decomposition. Returns a sample point
    in one (if return_all_samples=False) or all cells over which
    the system holds. If the system is invalid over all cells, then
    we return None.

    Parameters
    ==========

    seq: a list/tuple/set
        Listing all the (in)equalities that are needed to be solved
    gens: generators
        generators of the (in)equalities in seq for which we want the
        solutions
    strict: a boolean (default is False)
        if strict is True, NotImplementedError will be raised if
        the solution is known to be incomplete (which can occur if
        not all solutions are expressible in radicals)
    args: Keyword arguments
        Special options for solving the (in)equalities.


    Returns
    =======

    List[Tuple]
        a list of tuples with elements being solutions for the
        symbols in the order they were passed as gens
    None
        None is returned if no solutions exist.

    Examples
    ========

    ADD EXAMPLES HERE
    """
    # prepare the atoms
    atoms = [Poly(p.lhs - p.rhs, gens) for p in seq]

    sample_points = cylindrical_algebraic_decomposition(atoms, gens)
    valid_samples = []

    for sample in sample_points:
        if all([expr.subs(sample) for expr in seq]):
            valid_samples.append(sample)
            if return_one_sample:
                break
    
    return valid_samples






# HONG PROJECTOR OPERATOR AND OPERATIONS USED FOR IT

"""
Reducta

The reducta of a function f, treated as a univariate function of
main variable (mvar) x, is red(f) = f - lt(f), where lt is the
leading term.

The ith level reducta of f, red^i(f), is defined recursively.
red^0(f) = f
red^i(f) = red(red^{i-1}(f))

The reducta set RED(f) is defined as: {red^i(f) | 0 <= i <= deg(f)}
We exclude 0 from the set by definition as well.

"""

def red(f, mvar):
    return f - LT(f, mvar)

def red_level(f, mvar, level):
    if level == 0:
        return f
    else:
        red_curr = f
        for _ in range(level):
            red_curr = red(red_curr, mvar)
        return red_curr
    
def red_set(f, mvar):
    reds = []

    # handle this case here because otherwise sympy says deg= -infty
    if not f.free_symbols:
        return []
    
    for i in range(degree(f, mvar) + 1):
        if i == 0:
            red_curr = f
        else:
            red_curr = red(reds[i-1], mvar)
        
        if red_curr == 0:
            break
        else:
            reds.append(red_curr)
    
    return reds


"""
Principal subresultant coefficients

In Basu, Pollack, and Roy (2006) this is defined as a subdeterminant
of the Sylvester-Habicht matrix. However, here we can just use the
subresultants already built in to SymPy, and take coefficients.
"""

def psc_set(f, g, mvar):
    subres_polys = subresultants(f, g, mvar)[::-1]
    
    # subresultants returns it in decreasing order of degree
    # so we reverse it so that the ith subres is at most degree i
    # to match with the definition in Basu book
    # then the coefficient of ith subres of mvar^{i} is the ith psc

    final_psc = []

    for i in range(len(subres_polys)):
        p = subres_polys[i]
        if p == 0:
            continue
        
        curr_coeff = Poly(p, mvar).nth(i)
        if curr_coeff == 0:
            continue
        
        final_psc.append(curr_coeff)
    
    return final_psc


# gets real roots, numeric if not algebraic
def get_nice_roots(poly):
    # its a constant
    if not poly.free_symbols:
        return []
    
    try:
        roots = real_roots(poly)
    except NotImplementedError:
        return [root for root in nroots(poly) if root.is_real]

    
    roots_return = []
    for root in roots:
        """
        Want to avoid CRootOf (for now?)

        Note that we can't simply check for isinstance(ComplexRootOf)
        Sometimes, expressions can involve ComplexRootOf, eg:
        >>> Poly(-15*z**4 + 120*z**2 + 16).real_roots()[0]
        2*CRootOf(15*x**4 - 30*x**2 - 1, 0)

        To get around this, we use .has(). 
        """

        if root.has(ComplexRootOf):
            roots_return.append(root.evalf()) # not sure how precise?
        else:
            roots_return.append(root)

    return roots_return


# HONG PROJECTION OPERATOR (1990)
# input: set F of k-variate polynomials
# output: set F' of (k-1)-variate polynomials such that a CAD of R^{k-1} can be lifted to R^k

def projone(F, mvar):
    """
    Computes the PROJ1 operator as defined in Hong 1990.

    Let F be a set of polynomials with a given mvar. Then,

    PROJ1 = \cup_{f \in F, g \in RED(f)} ( ldcf(g) \cup PSC(g, D(g)))

    where RED is the reducta set, ldcf is the leading coefficient,
    PSC is the principal subresultant coefficient set, and D is the
    derivative operator.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.


    Returns
    =======

    Set
        A set of projection factors.

    Examples
    ========

    ADD EXAMPLES HERE
    """
    proj_set = set()
    for f in F:
        for g in red_set(f, mvar):
            proj_set.add(LC(g, mvar))
            proj_set.update(psc_set(g, diff(g, mvar), mvar))    
    
    return proj_set
            

def projtwo(F, mvar):
    """
    Computes the PROJ2* operator as defined in Hong (1990).
    This is an updated version of the PROJ2 operator from Collins.
    We will just call it PROJ2 here.

    Let F be a set of polynomials with a given mvar. Then,

    PROJ2 = \cup_{f,g \in F, f < g} \cup_{f' \in RED(f)} PSC(f', g)

    where RED is the reducta set, < indicates an arbitray "linear 
    ordering" to not loop over redundant pairs, and PSC is the 
    principal subresultant coefficient set.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    Set
        A set of projection factors.

    Examples
    ========

    ADD EXAMPLES HERE
    """

    proj_set = set()
    for i in range(len(F)):
        # impose "linear ordering"
        for j in range(i, len(F)):
            f, g = F[i], F[j]
            for f_ in red_set(f, mvar):
                proj_set.update(psc_set(f_, g, mvar))
    
    return proj_set


def hongproj(F, mvar):
    """
    The Hong projection operator, as defined in Hong(1990).
    PROJH takes a set of k-variate polynomials F, with an mvar, and  
    returns a set of (k-1)-variate polynomials F, with the mvar
    eliminated. These projection factors satisfy the property that a
    CAD of R^{k-1} can be lifted to a CAD of R^k.

    The Hong projector, PROJH, is defined as:

    PROJH(F) = PROJ1(F) \cup PROJ2(F)

    PROJ1 = \cup_{f \in F, g \in RED(f)} (ldcf(g) \cup PSC(g, D(g)))

    PROJ2 = \cup_{f,g \in F, f < g} \cup_{f' \in RED(f)} PSC(f', g)

    where RED is the reducta set, ldcf is the leading coefficient,
    PSC is the principal subresultant coefficient set, < indicates
    an arbitray "linear ordering" to not loop over redundant pairs,
    and D is the derivative operator.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    List
        The set of projection factors.
    """

    return list(projone(F, mvar).union(projtwo(F, mvar)))


def cylindrical_algebraic_decomposition(F, gens):
    """
    Calculates a cylindrical algebraic decomposition adapted to F.
    Uses the Hong projection operator. Returns sample points which
    represent cells over which each f \in F is sign-invariant. It
    projects iteratively down to lower-dimension spaces according to
    the list of generators given in gens.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    gens: a list of generators

    Returns
    =======

    ??: sample points
        
    """

    # Compute the projection sets
    projs_set = [F]
    for i in range(len(gens) - 1):
        projs_set.append(hongproj(projs_set[-1], gens[i]))

    # Lifting
    sample_points = [{}]

    for i in reversed(range(len(gens))):
        projs = projs_set[i]
        gen = gens[i]

        new_sample_points = []

        for point in sample_points:
            roots = set()
            for proj in projs:
                roots.update(get_nice_roots(proj.subs(point)))
            roots = sorted(list(roots), reverse=False)

            # Calculate sample points
            if not roots:
                samples = [0]
            elif len(roots) == 1:
                samples = [roots[0] - 1, roots[0], roots[0] + 1]
            else:
                samples = [roots[0] - 1]  # Point below the smallest root
                for r1, r2 in zip(roots, roots[1:]):
                    samples.extend([r1, (r1 + r2) / 2])
                samples.extend([roots[-1], roots[-1] + 1])  # Last root and point above it
            
            for value in samples:
                new_point = point.copy()
                new_point[gen] = value
                new_sample_points.append(new_point)
        
        sample_points = new_sample_points
    

    return sample_points
            


        








        

