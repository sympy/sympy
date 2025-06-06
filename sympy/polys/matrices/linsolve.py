#
# sympy.polys.matrices.linsolve module
#
# This module defines the _linsolve function which is the internal workhorse
# used by linsolve. This computes the solution of a system of linear equations
# using the SDM sparse matrix implementation in sympy.polys.matrices.sdm. This
# is a replacement for solve_lin_sys in sympy.polys.solvers which is
# inefficient for large sparse systems due to the use of a PolyRing with many
# generators:
#
#     https://github.com/sympy/sympy/issues/20857
#
# The implementation of _linsolve here handles:
#
# - Extracting the coefficients from the Expr/Eq input equations.
# - Constructing a domain and converting the coefficients to
#   that domain.
# - Using the SDM.rref, SDM.nullspace etc methods to generate the full
#   solution working with arithmetic only in the domain of the coefficients.
#
# The routines here are particularly designed to be efficient for large sparse
# systems of linear equations although as well as dense systems. It is
# possible that for some small dense systems solve_lin_sys which uses the
# dense matrix implementation DDM will be more efficient. With smaller systems
# though the bulk of the time is spent just preprocessing the inputs and the
# relative time spent in rref is too small to be noticeable.
#

from collections import defaultdict

from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.singleton import S

from sympy.polys.constructor import construct_domain
from sympy.polys.solvers import PolyNonlinearError

from .domainmatrix import DomainMatrix
from .exceptions import DMNonInvertibleMatrixError
from .sdm import (
    sdm_irref,
    sdm_particular_from_rref,
    sdm_nullspace_from_rref
)

from sympy.utilities.misc import filldedent


def _linsolve(eqs, syms):

    """Solve a linear system of equations.

    Examples
    ========

    Solve a linear system with a unique solution:

    >>> from sympy import symbols, Eq
    >>> from sympy.polys.matrices.linsolve import _linsolve
    >>> x, y = symbols('x, y')
    >>> eqs = [Eq(x + y, 1), Eq(x - y, 2)]
    >>> _linsolve(eqs, [x, y])
    {x: 3/2, y: -1/2}

    In the case of underdetermined systems the solution will be expressed in
    terms of the unknown symbols that are unconstrained:

    >>> _linsolve([Eq(x + y, 0)], [x, y])
    {x: -y, y: y}

    """
    # Convert to sparse augmented matrix (len(eqs) x (len(syms)+1))
    eqsdict, const = _linear_eq_to_dict(eqs, syms)
    Aaug = sympy_dict_to_dm(eqsdict, const, syms)
    return _linsolve_aug(Aaug, syms)


def _linsolve_aug(Aaug, syms):
    """Solve linear system represented as an augmented DomainMatrix.

    Examples
    ========

    >>> from sympy import symbols, Matrix, linear_eq_to_matrix
    >>> from sympy.polys.matrices.linsolve import _linsolve_aug
    >>> x, y = symbols('x, y')
    >>> eqs = [x + y - 1, x - y - 2]
    >>> A, b = linear_eq_to_matrix(eqs, [x, y])
    >>> Aaug = Matrix.hstack(A, b).to_DM()
    >>> Aaug
    DomainMatrix({0: {0: 1, 1: 1, 2: 1}, 1: {0: 1, 1: -1, 2: 2}}, (2, 3), ZZ)
    >>> _linsolve_aug(Aaug, [x, y])
    {x: 3/2, y: -1/2}
    """
    K = Aaug.domain

    try:
        P, den, V, free_variables = Aaug[:,:-1].solve_den_general(Aaug[:,-1])
    except DMNonInvertibleMatrixError:
        return None

    Kf = K

    if not K.is_one(den):
        if not K.is_Field:
            Kf = K.get_field()
            P = P.convert_to(Kf)
            V = V.convert_to(Kf)
        P = P / den
        V = V / den

    P = P.transpose().to_dod().get(0, {})
    Vd = V.to_dod()
    assert len(Vd) == len(free_variables)
    V = [Vd[i] for i in range(len(free_variables))]

    # No solution:
    if P is None:
        return None

    # Collect together terms from particular and nullspace:
    sol = defaultdict(list)
    for i, v in P.items():
        sol[syms[i]].append(Kf.to_sympy(v))
    for npi, Vi in zip(free_variables, V):
        sym = syms[npi]
        for i, v in Vi.items():
            sol[syms[i]].append(sym * Kf.to_sympy(v))

    # Use a single call to Add for each term:
    sol = {s: Add(*terms) for s, terms in sol.items()}

    # Fill in the zeros:
    zero = S.Zero
    for s in set(syms) - set(sol):
        sol[s] = zero

    # All done!
    return sol


def _particular_nullspace(Aaug):
    """Return a particular solution and nullspace basis for the
    augmented matrix ``Aaug``. The augmented matrix is assumed to
    be in reduced row echelon form. The particular solution is
    returned as a dictionary mapping column indices to the
    corresponding values.
    """
    # Number of unknowns (columns in the non-augmented matrix)
    nsyms = Aaug.shape[1] - 1

    K = Aaug.domain

    # sdm_irref has issues with float matrices. This uses the ddm_rref()
    # function. When sdm_rref() can handle float matrices reasonably this
    # should be removed...
    if K.is_RealField or K.is_ComplexField:
        Aaug = Aaug.to_ddm().rref()[0].to_sdm()

    # Compute reduced-row echelon form (RREF)
    Arref, pivots, nzcols = sdm_irref(Aaug)

    # No solution:
    if pivots and pivots[-1] == nsyms:
        return None, None, None

    # Particular solution for non-homogeneous system:
    P = sdm_particular_from_rref(Arref, nsyms+1, pivots)

    # Nullspace - general solution to homogeneous system
    # Note: using nsyms not nsyms+1 to ignore last column
    V, free_variables = sdm_nullspace_from_rref(Arref, K.one, nsyms, pivots, nzcols)

    return P, V, free_variables


def sympy_dict_to_dm(eqs_coeffs, eqs_rhs, syms):
    """Convert a system of dict equations to a sparse augmented matrix"""
    elems = set(eqs_rhs).union(*(e.values() for e in eqs_coeffs))
    K, elems_K = construct_domain(elems, field=True, extension=True)
    elem_map = dict(zip(elems, elems_K))

    neqs = len(eqs_coeffs)
    nsyms = len(syms)
    sym2index = dict(zip(syms, range(nsyms)))
    eqsdict = []
    for eq, rhs in zip(eqs_coeffs, eqs_rhs):
        eqdict = {sym2index[s]: elem_map[c] for s, c in eq.items()}
        if rhs:
            eqdict[nsyms] = -elem_map[rhs]
        if eqdict:
            eqsdict.append(eqdict)

    dod = dict(enumerate(eqsdict))
    shape = (neqs, nsyms + 1)
    Aaug = DomainMatrix.from_dod(dod, shape, K)

    return Aaug


def _linear_eq_to_dict(eqs, syms):
    """Convert a system Expr/Eq equations into dict form, returning
    the coefficient dictionaries and a list of syms-independent terms
    from each expression in ``eqs```.

    Examples
    ========

    >>> from sympy.polys.matrices.linsolve import _linear_eq_to_dict
    >>> from sympy.abc import x
    >>> _linear_eq_to_dict([2*x + 3], {x})
    ([{x: 2}], [3])
    """
    coeffs = []
    ind = []
    symset = set(syms)
    for e in eqs:
        if e.is_Equality:
            coeff, terms = _lin_eq2dict(e.lhs, symset)
            cR, tR = _lin_eq2dict(e.rhs, symset)
            # there were no nonlinear errors so now
            # cancellation is allowed
            coeff -= cR
            for k, v in tR.items():
                if k in terms:
                    terms[k] -= v
                else:
                    terms[k] = -v
            # don't store coefficients of 0, however
            terms = {k: v for k, v in terms.items() if v}
            c, d = coeff, terms
        else:
            c, d = _lin_eq2dict(e, symset)
        coeffs.append(d)
        ind.append(c)
    return coeffs, ind


def _lin_eq2dict(a, symset):
    """return (c, d) where c is the sym-independent part of ``a`` and
    ``d`` is an efficiently calculated dictionary mapping symbols to
    their coefficients. A PolyNonlinearError is raised if non-linearity
    is detected.

    The values in the dictionary will be non-zero.

    Examples
    ========

    >>> from sympy.polys.matrices.linsolve import _lin_eq2dict
    >>> from sympy.abc import x, y
    >>> _lin_eq2dict(x + 2*y + 3, {x, y})
    (3, {x: 1, y: 2})
    """
    if a in symset:
        return S.Zero, {a: S.One}
    elif a.is_Add:
        terms_list = defaultdict(list)
        coeff_list = []
        for ai in a.args:
            ci, ti = _lin_eq2dict(ai, symset)
            coeff_list.append(ci)
            for mij, cij in ti.items():
                terms_list[mij].append(cij)
        coeff = Add(*coeff_list)
        terms = {sym: Add(*coeffs) for sym, coeffs in terms_list.items()}
        return coeff, terms
    elif a.is_Mul:
        terms = terms_coeff = None
        coeff_list = []
        for ai in a.args:
            ci, ti = _lin_eq2dict(ai, symset)
            if not ti:
                coeff_list.append(ci)
            elif terms is None:
                terms = ti
                terms_coeff = ci
            else:
                # since ti is not null and we already have
                # a term, this is a cross term
                raise PolyNonlinearError(filldedent('''
                    nonlinear cross-term: %s''' % a))
        coeff = Mul._from_args(coeff_list)
        if terms is None:
            return coeff, {}
        else:
            terms = {sym: coeff * c for sym, c in terms.items()}
            return  coeff * terms_coeff, terms
    elif not a.has_xfree(symset):
        return a, {}
    else:
        raise PolyNonlinearError('nonlinear term: %s' % a)
