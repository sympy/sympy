"""Low-level linear systems solver. """

from __future__ import print_function, division

from sympy.utilities.iterables import connected_components

from sympy.matrices import MutableDenseMatrix
from sympy.polys.domains import EX
from sympy.polys.rings import sring
from sympy.polys.polyerrors import NotInvertible
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

def eqs_to_matrix2(eqs_coeffs, eqs_rhs, gens, domain):
    sym2index = {x: n for n, x in enumerate(gens)}
    nrows = len(eqs_coeffs)
    ncols = len(gens) + 1
    rows = [[domain.zero] * ncols for _ in range(nrows)]
    for row, eq_coeff, eq_rhs in zip(rows, eqs_coeffs, eqs_rhs):
        for sym, coeff in eq_coeff.items():
            row[sym2index[sym]] = coeff
        row[-1] = -eq_rhs

    M = DomainMatrix(rows, (nrows, ncols), domain)

    return M

def sympy_eqs_to_ring(eqs, symbols):
    try:
        K, eqs_K = sring(eqs, symbols, field=True, extension=True)
    except NotInvertible:
        # https://github.com/sympy/sympy/issues/18874
        K, eqs_K = sring(eqs, symbols, domain=EX)
    return eqs_K, K.to_domain()

def solve_lin_sys(eqs, ring, _raw=True):
    """Solve a system of linear equations.

    If ``_raw`` is False, the keys and values in the returned dictionary
    will be of type Expr (and the unit of the field will be removed from
    the keys) otherwise the low-level polys types will be returned, e.g.
    PolyElement: PythonRational.
    """
    as_expr = not _raw

    assert ring.domain.is_Field

    eqs_dict = [dict(eq) for eq in eqs]

    one_monom = ring.one.monoms()[0]
    zero = ring.domain.zero

    eqs_rhs = []
    eqs_coeffs = []
    for eq_dict in eqs_dict:
        eq_rhs = eq_dict.pop(one_monom, zero)
        eq_coeffs = {ring.gens[monom.index(1)]:coeff for monom, coeff in eq_dict.items()}
        if not eq_coeffs:
            if not eq_rhs:
                continue
            else:
                return None
        eqs_rhs.append(eq_rhs)
        eqs_coeffs.append(eq_coeffs)

    result = _solve_lin_sys(eqs_coeffs, eqs_rhs, ring)

    if result is not None and as_expr:

        def to_sympy(x):
            as_expr = getattr(x, 'as_expr', None)
            if as_expr:
                return as_expr()
            else:
                return ring.domain.to_sympy(x)

        tresult = {to_sympy(sym): to_sympy(val) for sym, val in result.items()}

        # Remove 1.0x
        result = {}
        for k, v in tresult.items():
            if k.is_Mul:
                c, s = k.as_coeff_Mul()
                result[s] = v/c
            else:
                result[k] = v

    return result

def _solve_lin_sys(eqs_coeffs, eqs_rhs, ring):

    V = ring.gens
    E = []
    for eq_coeffs in eqs_coeffs:
        syms = list(eq_coeffs)
        E.extend(zip(syms[:-1], syms[1:]))
    G = V, E

    components = connected_components(G)

    sym2comp = {}
    for n, component in enumerate(components):
        for sym in component:
            sym2comp[sym] = n

    subsystems = [([], []) for _ in range(len(components))]
    for eq_coeff, eq_rhs in zip(eqs_coeffs, eqs_rhs):
        sym = next(iter(eq_coeff), None)
        sub_coeff, sub_rhs = subsystems[sym2comp[sym]]
        sub_coeff.append(eq_coeff)
        sub_rhs.append(eq_rhs)

    sol = {}
    for subsystem in subsystems:
        subsol = _solve_lin_sys_component(subsystem[0], subsystem[1], ring)
        if subsol is None:
            return None
        sol.update(subsol)

    return sol


def _solve_lin_sys_component(eqs_coeffs, eqs_rhs, ring):

    # transform from equations to matrix form
    matrix = eqs_to_matrix2(eqs_coeffs, eqs_rhs, ring.gens, ring.domain)

    # solve by row-reduction
    echelon, pivots = matrix.rref()

    # construct the returnable form of the solutions
    keys = ring.gens

    if pivots and pivots[-1] == len(keys):
        return None

    if len(pivots) == len(keys):
        sol = []
        for s in [row[-1] for row in echelon.rows]:
            a = s
            sol.append(a)
        sols = dict(zip(keys, sol))
    else:
        sols = {}
        g = ring.gens
        echelon = echelon.rows
        for i, p in enumerate(pivots):
            v = echelon[i][-1] - sum(echelon[i][j]*g[j] for j in range(p+1, len(g)))
            v = v
            sols[keys[p]] = v

    return sols
