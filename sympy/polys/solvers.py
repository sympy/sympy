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

def solve_lin_sys(eqs, ring, iszerofunc=lambda x: not x, scalefunc=None, elimfunc=None):
    """
    Solve a system of linear equations over a field.

    The optional arguments ``iszerofunc``, ``scalefunc`` and ``elimfunc``
    can to be used to specify dfferent operations than the ones given by
    ``ring``. ``iszerofunc`` will be used to decide if a domain element is
    zero, ``scalefunc`` will be used to scale the row of the current pivot
    element and ``elimfunc`` will be used in the elimination step.

    For example, a linear system of integer equations can be solved over
    `\mathbb Z_p` by defining those functions as:

        |  iszerofunc = lambda x: not (x % p)
        |  scalefunc = lambda x, _, scale: (x * invert(scale, p)) % p
        |  elimfunc = lambda x, y, scale: (x - scale*y) % p

    Note that ``scalefunc`` and ``elimfunc`` have to take three arguments,
    where the third one has to be named ``scale``.

    """
    # transform from equations to matrix form
    matrix = eqs_to_matrix(eqs, ring)

    # solve by row-reduction
    echelon, pivots = matrix.rref(iszerofunc=iszerofunc, scalefunc=scalefunc, elimfunc=elimfunc)

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
