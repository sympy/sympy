
from sympy.polys import Poly, LexPoly, roots, \
    SymbolsError, PolynomialError
from sympy.polys.algorithms import poly_groebner

from sympy.utilities import any, all

def solve_poly_system(system, *symbols):
    """Solves a system of polynomial equations.

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

       >>> from sympy import *
       >>> from sympy.abc import x, y

       >>> solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y)
       [(0, 0), (2, -2**(1/2)), (2, 2**(1/2))]

       For more information on the implemented algorithm refer to:

       [1] B. Buchberger, Groebner Bases: A Short Introduction for
           Systems Theorists,  In: R. Moreno-Diaz,  B. Buchberger,
           J.L. Freire, Proceedings of EUROCAST'01, February, 2001

       [2] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 112

    """
    def is_univariate(f):
        """Returns True if 'f' is univariate in its last variable. """
        for monom in f.iter_monoms():
            if any(exp > 0 for exp in monom[:-1]):
                return False

        return True

    def solve_reduced_system(system, entry=False):
        """Recursively solves reduced polynomial systems. """
        basis = poly_groebner(system)

        if len(basis) == 1 and basis[0].is_one:
            if not entry:
                return []
            else:
                return None

        univariate = filter(is_univariate, basis)

        if len(univariate) == 1:
            f = univariate.pop()
        else:
            raise PolynomialError("Not a zero-dimensional system")

        zeros = roots(Poly(f, f.symbols[-1])).keys()

        if not zeros:
            return []

        if len(basis) == 1:
            return [ [zero] for zero in zeros ]

        solutions = []

        for zero in zeros:
            new_system = []

            for poly in basis[:-1]:
                eq = poly.evaluate((poly.symbols[-1], zero))

                if not eq.is_zero:
                    new_system.append(eq)

            for solution in solve_reduced_system(new_system):
                solutions.append(solution + [zero])

        return solutions

    if hasattr(system, "__iter__"):
        system = list(system)
    else:
        raise TypeError("Expected iterable container, got %s" % system)

    f = system.pop(0)

    if not isinstance(f, Poly):
        f = LexPoly(f, *symbols)
    else:
        if not symbols:
            f = LexPoly(f)
        else:
            raise SymbolsError("Redundant symbols were given")

    head, tail = f.unify_with(system)

    solutions = solve_reduced_system([head] + tail, True)

    if solutions is None:
        return None
    else:
        return sorted(tuple(s) for s in solutions)

