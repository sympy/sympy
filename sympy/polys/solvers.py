"""Low-level linear systems solver. """

from sympy.matrices import Matrix, zeros

class RawMatrix(Matrix):
    _sympify = staticmethod(lambda x: x)

def solve_lin_sys(eqs, ring):
    """Solve a system of linear equations. """

    # transform from equations to matrix form
    xs = ring.gens
    M = zeros(len(eqs), len(xs)+1, cls=RawMatrix)

    for j, e_j in enumerate(eqs):
        for i, x_i in enumerate(xs):
            M[j, i] = e_j.coeff(x_i)
        M[j, -1] = -e_j.coeff(1)

    eqs = M

    # solve by row-reduction
    eschelon, pivots = eqs.rref(iszerofunc=lambda x: not x, simplify=lambda x: x)

    # construct the returnable form of the solutions
    p = len(pivots)
    if p > len(xs):
        return None

    sols = eschelon[:p, p:]*RawMatrix([ [x] for x in xs[p:] ] + [[1]])
    return dict(zip(xs, sols))
