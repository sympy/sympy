"""Low-level linear systems solver. """

from __future__ import print_function, division

from sympy.matrices import MutableDenseMatrix
from sympy.polys.polymatrix import DomainMatrix

class RawMatrix(MutableDenseMatrix):
    _sympify = staticmethod(lambda x: x)

def eqs_to_matrix(eqs, ring):
    """Transform from equations to matrix form. """
    xs = ring.gens

    shape = (len(eqs), len(xs) + 1)

    rows = [[e_j.coeff(x_i) for x_i in xs] for e_j in eqs]
    for row_j, e_j in zip(rows, eqs):
        row_j.append(-e_j.coeff(1))

    M = DomainMatrix(rows, shape, ring.domain)

    return M

def solve_lin_sys(eqs, ring, _raw=True):
    """Solve a system of linear equations.

    If ``_raw`` is False, the keys and values in the returned dictionary
    will be of type Expr (and the unit of the field will be removed from
    the keys) otherwise the low-level polys types will be returned, e.g.
    PolyElement: PythonRational.
    """
    as_expr = not _raw

    assert ring.domain.is_Field

    # transform from equations to matrix form
    matrix = eqs_to_matrix(eqs, ring)

    # solve by row-reduction
    echelon, pivots = matrix.rref()

    # construct the returnable form of the solutions
    keys = ring.symbols if as_expr else ring.gens

    if pivots[-1] == len(keys):
        return None

    if len(pivots) == len(keys):
        sol = []
        for s in [row[-1] for row in echelon.rows]:
            a = ring.ground_new(s)
            if as_expr:
                a = a.as_expr()
            sol.append(a)
        sols = dict(zip(keys, sol))
    else:
        sols = {}
        g = ring.gens
        echelon = echelon.rows
        for i, p in enumerate(pivots):
            v = echelon[i][-1] - sum(echelon[i][j]*g[j] for j in range(p+1, len(g)))
            v = (v + ring.zero)
            if as_expr:
                v = v.as_expr()
            sols[keys[p]] = v

    return sols
