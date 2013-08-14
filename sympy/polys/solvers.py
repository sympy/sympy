"""Low-level linear systems solver. """

from __future__ import print_function, division

from sympy.matrices import Matrix, zeros

class RawMatrix(Matrix):
    _sympify = staticmethod(lambda x: x)

def eqs_to_matrix(eqs, ring):
    """Transform from equations to matrix form. """
    xs = ring.gens
    M = zeros(len(eqs), len(xs)+1, cls=RawMatrix)

    for j, e_j in enumerate(eqs):
        for i, x_i in enumerate(xs):
            M[j, i] = e_j.coeff(x_i)
        M[j, -1] = -e_j.coeff(1)

    return M

def solve_lin_sys(eqs, ring):
    """Solve a system of linear equations. """
    assert ring.domain.has_Field

    # transform from equations to matrix form
    matrix = eqs_to_matrix(eqs, ring)

    # solve by row-reduction
    echelon, pivots = matrix.rref(iszerofunc=lambda x: not x, simplify=lambda x: x)

    # construct the returnable form of the solutions
    xs = ring.gens

    if pivots[-1] == len(xs):
        return None
    elif len(pivots) == len(xs):
        sol = [ ring.ground_new(s) for s in echelon[:, -1] ]
        return dict(zip(xs, sol))
    else:
        sols = {}
        for i, p in enumerate(pivots):
            vect = RawMatrix([ [-x] for x in xs[p+1:] ] + [[ring.one]])
            sols[xs[p]] = (echelon[i, p+1:]*vect)[0]

        return sols
