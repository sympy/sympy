# Algorithms for computing the reduced row echelon form of a matrix.
#
# We need to choose carefully which algorithms to use depending on the domain,
# shape, and sparsity of the matrix as well as things like the bit count in the
# case of ZZ or QQ. This is important because the algorithms have different
# performance characteristics in the extremes of dense vs sparse.
#
# In all cases we use the sparse implementations but we need to choose between
# Gauss-Jordan elimination with division and fraction-free Gauss-Jordan
# elimination. For very sparse matrices over ZZ with low bit counts it is
# asymptotically faster to use Gauss-Jordan elimination with division. For
# dense matrices with high bit counts it is asymptotically faster to use
# fraction-free Gauss-Jordan.
#
# The most important thing is to get the extreme cases right because it can
# make a big difference. In between the extremes though we have to make a
# choice and here we use empirically determined thresholds based on timings
# with random sparse matrices.
#
# In the case of QQ we have to consider the denominators as well. If the
# denominators are small then it is faster to clear them and use fraction-free
# Gauss-Jordan over ZZ. If the denominators are large then it is faster to use
# Gauss-Jordan elimination with division over QQ.

from sympy.polys.matrices.exceptions import DMNotAField

from sympy.polys.domains import ZZ, QQ

from sympy.polys.matrices.sdm import SDM, sdm_irref, sdm_rref_den
from sympy.polys.matrices.ddm import DDM
from sympy.polys.matrices.dense import ddm_irref


def dm_rref(M):
    """
    Compute the reduced row echelon form of a ``DomainMatrix``.

    Choose the best algorithm depending on the domain, shape, and sparsity of
    the matrix as well as things like the bit count in the case of ZZ or QQ.

    Examples
    ========

    >>> from sympy import QQ
    >>> from sympy.polys.matrices import DM
    >>> from sympy.polys.matrices.rref import dm_rref
    >>> M = DM([[1, 2], [3, 4]], QQ)
    >>> M_rref, pivots = dm_rref(M)
    >>> M_rref.to_Matrix()
    Matrix([
    [1, 0],
    [0, 1]])
    >>> pivots
    (0, 1)

    See Also
    ========

    sympy.polys.matrices.domainmatrix.DomainMatrix.rref
        The ``DomainMatrix`` method that calls this function.
    sympy.polys.matrices.rref.dm_rref_den
        Alternative function for computing RREF with denominator.
    """
    K = M.domain

    if not K.is_Field:
        raise DMNotAField("Use .rref_den() or convert to a field with .to_field()")

    dense = M.rep.fmt == 'dense'

    # Do not switch to the sparse implementation for EX because the domain does
    # not have proper canonicalization and the sparse implementation gives
    # equivalent but non-identical results over EX from performing arithmetic
    # in a different order. Specifically test_issue_23718 ends up getting a
    # more complicated expression when using the sparse implementation.
    # Probably the best fix for this is something else but for now we stick
    # with the dense implementation for EX if the matrix is already dense.
    if K.is_EX and dense:
        return _dm_rref_gj_div_dense(M)

    # The sparse implementations are always faster when we can use them.
    if dense:
        M = M.to_sparse()

    if K.is_QQ:
        M_rref, pivots = _dm_rref_QQ(M)
    else:
        # This is definitely not optimal. More work needs to be done to
        # determine the best algorithm for each domain. Generally it is better
        # to use rref_den for nontrivial field domains in any case though.
        M_rref, pivots = _dm_rref_gj_div(M)

    if dense:
        M_rref = M_rref.to_dense()

    return M_rref, pivots


def dm_rref_den(M, keep_domain=True):
    """
    Compute the reduced row echelon form of a matrix with denominator.

    Examples
    ========

    >>> from sympy import ZZ
    >>> from sympy.polys.matrices import DM
    >>> from sympy.polys.matrices.rref import dm_rref_den
    >>> M = DM([[1, 2], [3, 4]], ZZ)
    >>> M_rref, den, pivots = dm_rref_den(M)
    >>> M_rref.to_Matrix()
    Matrix([
    [-2,  0],
    [ 0, -2]])
    >>> den
    -2
    >>> pivots
    (0, 1)

    Parameters
    ==========

    keep_domain : bool, optional
        If ``True`` then the domain of the matrix is preserved. If ``False``
        then the domain might be preserved or it might be changed to an
        associated field or ring (if that is potentially faster).

    See Also
    ========

    sympy.polys.matrices.domainmatrix.DomainMatrix.rref_den
        The ``DomainMatrix`` method that calls this function.
    """
    K = M.domain

    # The sparse implementations are always faster
    dense = M.rep.fmt == 'dense'
    if dense:
        M = M.to_sparse()

    if K.is_ZZ:
        M_rref, den, pivots = _dm_rref_den_ZZ(M, keep_domain=keep_domain)
    elif K.is_QQ:
        M_rref, den, pivots = _dm_rref_den_QQ(M, keep_domain=keep_domain)
    else:
        M_rref, den, pivots = _dm_rref_den_gj_ff(M)

    if dense:
        M_rref = M_rref.to_dense()

    return M_rref, den, pivots


def _dm_rref_gj_div(M):
    """Compute RREF using sparse Gauss-Jordan elimination with division."""
    M_rref_d, pivots, _ = sdm_irref(M.rep)
    M_rref_sdm = SDM(M_rref_d, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_sdm), pivots


def _dm_rref_gj_div_dense(M):
    """Compute RREF using dense Gauss-Jordan elimination with division."""
    ddm = M.rep.copy()
    pivots = ddm_irref(ddm)
    M_rref_ddm = DDM(ddm, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_ddm), pivots


def _dm_rref_den_gj_ff(M):
    """Compute RREF using sparse fraction-free Gauss-Jordan elimination."""
    M_rref_d, den, pivots = sdm_rref_den(M.rep, M.domain)
    M_rref_sdm = SDM(M_rref_d, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_sdm), den, pivots


def _dm_rref_den_ZZ(Mz, keep_domain=True):
    """Compute RREF over ZZ using the fastest method."""
    # Choose the fastest method
    method = _dm_rref_den_ZZ_fastest_method(Mz)

    if method == 'rref_QQ':
        # Use Gauss-Jordan over QQ for sparse matrices.
        M_rref, pivots = _dm_rref_gj_div(Mz.convert_to(QQ))

        if keep_domain:
            _, M_rref = M_rref.clear_denoms(convert=True)

        if pivots:
            den = M_rref[0, pivots[0]].element
        else:
            den = ZZ.one

    elif method == 'rref_den_ZZ':
        # Otherwise use fraction-free Gauss-Jordan over ZZ.
        M_rref, den, pivots = _dm_rref_den_gj_ff(Mz)

    else:
        assert False  # pragma: no cover

    return M_rref, den, pivots


def _dm_rref_QQ(Mq):
    """Compute RREF over QQ using the fastest method."""
    # Choose the fastest method
    method = _dm_rref_QQ_fastest_method(Mq)

    if method == 'rref_QQ':
        # Sparse or has large denominators, use Gauss-Jordan over QQ.
        M_rref_q, pivots = _dm_rref_gj_div(Mq)

    elif method == 'rref_ZZ_clear_denoms':
        # Clear small denominators and use fraction-free Gauss-Jordan over ZZ.
        _, Mz = Mq.clear_denoms(convert=True)
        M_rref_z, den, pivots = _dm_rref_den_gj_ff(Mz)
        M_rref_q = M_rref_z.convert_to(QQ) / den

    elif method == 'rref_den_QQ':
        # Use fraction-free Gauss-Jordan over QQ
        M_rref_q, den, pivots = _dm_rref_den_gj_ff(Mq)
        M_rref_q /= den

    else:
        assert False  # pragma: no cover

    return M_rref_q, pivots


def _dm_rref_den_QQ(Mq, keep_domain=True):
    """Compute RREF over QQ using the fastest method."""
    # Choose the fastest method
    method = _dm_rref_QQ_fastest_method(Mq)

    if method == 'rref_QQ':
        # Use Gauss-Jordan with division over QQ.
        M_rref, pivots = _dm_rref_gj_div(Mq)
        den = QQ.one

    elif method == 'rref_den_QQ':
        # Use fraction-free Gauss-Jordan over QQ.
        M_rref, den, pivots = _dm_rref_den_gj_ff(Mq)

    elif method == 'rref_ZZ_clear_denoms':
        # Clear small denominators and use fraction-free Gauss-Jordan over ZZ.
        _, Mz = Mq.clear_denoms(convert=True)
        M_rref, den, pivots = _dm_rref_den_gj_ff(Mz)

        if keep_domain:
            M_rref = M_rref.convert_to(QQ)
            den = QQ.convert_from(den, ZZ)

    else:
        assert False  # pragma: no cover

    return M_rref, den, pivots


def _dm_rref_den_ZZ_fastest_method(Mz):
    """Return True if rref_QQ should be used."""
    # In the extreme of very sparse matrices and low bit counts it is faster to
    # use Gauss-Jordan elimination over QQ rather than fraction-free
    # Gauss-Jordan over ZZ. In the opposite extreme of dense matrices and high
    # bit counts it is faster to use fraction-free Gauss-Jordan over ZZ. These
    # two extreme cases need to be handled differently because they lead to
    # different asymptotic complexities. In between these two extremes we need
    # a threshold for deciding which method to use. This threshold is
    # determined empirically by timing the two methods with random matrices.

    # The disadvantage of using empirical timings is that future optimisations
    # might change the relative speeds so this can easily become out of date.
    # The main thing is to get the asymptotic complexity right for the extreme
    # cases though so the precise value of the threshold is hopefully not too
    # important.

    # Empirically determined parameter.
    PARAM = 10000

    # First compute the density. This is the average number of non-zero entries
    # per row but only counting rows that have at least one non-zero entry
    # since RREF can ignore fully zero rows.
    density, nrows_nz, ncols = _dm_row_density(Mz)

    # For small matrices use QQ if more than half the entries are zero.
    if nrows_nz < 10:
        if density < ncols/2:
            return 'rref_QQ'
        else:
            return 'rref_den_ZZ'

    # These are just shortcuts for the formula below.
    if density < 5:
        return 'rref_QQ'
    elif density > 5 + PARAM/nrows_nz:
        return 'rref_den_ZZ'  # pragma: no cover

    # Maximum bitsize of any entry.
    elements = _dm_elements(Mz)
    bits = max([e.bit_length() for e in elements], default=1)

    # Wideness parameter. This is 1 for square or tall matrices but >1 for wide
    # matrices.
    wideness = max(1, 2/3*ncols/nrows_nz)

    max_density = (5 + PARAM/(nrows_nz*bits**2)) * wideness

    if density < max_density:
        return 'rref_QQ'
    else:
        return 'rref_den_ZZ'


def _dm_rref_QQ_fastest_method(Mq):
    """Choose the fastest method for computing RREF over QQ."""
    # The same sorts of considerations apply here as in the case of ZZ. Here
    # though a new more significant consideration is what sort of denominators
    # we have and what to do with them so we focus on that.

    # First compute the density. This is the average number of non-zero entries
    # per row but only counting rows that have at least one non-zero entry
    # since RREF can ignore fully zero rows.
    density, _, ncols = _dm_row_density(Mq)

    # For sparse matrices use Gauss-Jordan elimination over QQ regardless.
    if density < min(5, ncols/2):
        return 'rref_QQ'

    # Compare the bit-length of the lcm of the denominators to the bit length
    # of the numerators.
    #
    # The threshold here is empirical: we prefer rref over QQ if clearing
    # denominators would result in a numerator matrix having 5x the bit size of
    # the current numerators.
    numers, denoms = _dm_QQ_numers_denoms(Mq)
    numer_bits = max([n.bit_length() for n in numers], default=1)

    denom_lcm = ZZ.one
    for d in denoms:
        denom_lcm = ZZ.lcm(denom_lcm, d)
        if denom_lcm.bit_length() > 5*numer_bits:
            return 'rref_QQ'

    # If we get here then the matrix is dense and the lcm of the denominators
    # is not too large compared to the numerators. For particularly small
    # denominators it is fastest just to clear them and use fraction-free
    # Gauss-Jordan over ZZ. With very small denominators this is a little
    # faster than using rref_den over QQ but there is an intermediate regime
    # where rref_den over QQ is significantly faster. The small denominator
    # case is probably very common because small fractions like 1/2 or 1/3 are
    # often seen in user inputs.

    if denom_lcm.bit_length() < 50:
        return 'rref_ZZ_clear_denoms'
    else:
        return 'rref_den_QQ'


def _dm_row_density(M):
    """Density measure for sparse matrices.

    Defines the "density", ``d`` as the average number of non-zero entries per
    row except ignoring rows that are fully zero. RREF can ignore fully zero
    rows so they are excluded. By definition ``d >= 1`` except that we define
    ``d = 0`` for the zero matrix.

    Returns ``(density, nrows_nz, ncols)`` where ``nrows_nz`` counts the number
    of nonzero rows and ``ncols`` is the number of columns.
    """
    # Assumes the SDM dict-of-dicts representation.
    ncols = M.shape[1]
    rows_nz = M.rep.values()
    if not rows_nz:
        return 0, 0, ncols
    else:
        nrows_nz = len(rows_nz)
        density = sum(map(len, rows_nz)) / nrows_nz
        return density, nrows_nz, ncols


def _dm_elements(M):
    """Return nonzero elements of a DomainMatrix."""
    elements, _ = M.to_flat_nz()
    return elements


def _dm_QQ_numers_denoms(Mq):
    """Returns the numerators and denominators of a DomainMatrix over QQ."""
    elements = _dm_elements(Mq)
    numers = [e.numerator for e in elements]
    denoms = [e.denominator for e in elements]
    return numers, denoms
