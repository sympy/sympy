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
#
# Timings for the various algorithms can be found at
#
#   https://github.com/sympy/sympy/issues/25410
#   https://github.com/sympy/sympy/pull/25443

from sympy.polys.domains import ZZ

from sympy.polys.matrices.sdm import SDM, sdm_irref, sdm_rref_den
from sympy.polys.matrices.ddm import DDM
from sympy.polys.matrices.dense import ddm_irref, ddm_irref_den


def dm_rref(M, *, method='auto'):
    """
    Compute the reduced row echelon form of a ``DomainMatrix``.

    Chooses the best algorithm depending on the domain, shape, and sparsity of
    the matrix as well as things like the bit count in the case of ZZ or QQ.
    The result is returned over the field associated with the domain of the
    matrix.

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

    Parameters
    ==========

    method : str, optional (``'auto'``, ``'GJ'``, ``'FF'``, ``'CD'``,
        ``'GJ_dense'``, ``'FF_dense'``, ``'CD_dense'``)

        The algorithm to use. The default is ``'auto'`` which chooses the best
        algorithm depending on the domain, shape, and sparsity of the matrix as
        well as properties of the elements e.g. the bit count in the case of
        :ref:`ZZ` or :ref:`QQ`.

        With ``method='GJ'`` the matrix is converted to the associated field
        domain and Gauss-Jordan elimination with division is used. With
        ``method='FF'`` fraction-free Gauss-Jordan elimination is used in the
        current domain and the result is converted to the associated field
        domain at the end. With ``method='CD'`` the denominators are cleared
        and fraction-free Gauss-Jordan elimination is used in the associated
        ring before converting to the associated field domain and dividing at
        the end.

        By default the sparse implementations of the different algorithms are
        used. To use the dense implementations instead use
        ``method='GJ_dense'``, ``method='FF_dense'``, or ``method='CD_dense'``.
        With ``method='auto'`` only the sparse implementations are used.

        In all cases the result is returned over the field associated with the
        domain of the matrix and the format of the returned matrix (sparse or
        dense) is always the same as the format of the input matrix.

    Returns
    =======

    M_rref : DomainMatrix
        The reduced row echelon form of the matrix over the field associated
        with the domain of the matrix.
    pivots : tuple
        The indices of the pivot columns.

    See Also
    ========

    sympy.polys.matrices.domainmatrix.DomainMatrix.rref
        The ``DomainMatrix`` method that calls this function.
    sympy.polys.matrices.rref.dm_rref_den
        Alternative function for computing RREF with denominator.
    """
    K = M.domain

    is_dense = M.rep.fmt == 'dense'

    # Do not switch to the sparse implementation for EX because the domain does
    # not have proper canonicalization and the sparse implementation gives
    # equivalent but non-identical results over EX from performing arithmetic
    # in a different order. Specifically test_issue_23718 ends up getting a
    # more complicated expression when using the sparse implementation.
    # Probably the best fix for this is something else but for now we stick
    # with the dense implementation for EX if the matrix is already dense.
    if method == 'auto' and K.is_EX and is_dense:
        return _dm_rref_GJ_dense(M)

    if method == 'auto':
        method, use_dense = _dm_rref_choose_method(M)
    else:
        use_dense = method.endswith('_dense')
        if use_dense:
            method = method[:-len('_dense')]

    if use_dense:
        if not is_dense:
            M = M.to_dense()

        M_rref, pivots = _dm_rref(M, method, use_dense)

        if not is_dense:
            M_rref = M_rref.to_sparse()

    else:
        if is_dense:
            M = M.to_sparse()

        M_rref, pivots = _dm_rref(M, method, use_dense)

        if is_dense:
            M_rref = M_rref.to_dense()

    return M_rref, pivots


def dm_rref_den(M, *, method='auto', keep_domain=True):
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

    method : str, optional (``'auto'``, ``'GJ'``, ``'FF'``, ``'CD'``,
        ``'GJ_dense'``, ``'FF_dense'``, ``'CD_dense'``)

        The algorithm to use. The default is ``'auto'`` which chooses the best
        algorithm depending on the domain, shape, and sparsity of the matrix as
        well as properties of the elements e.g. the bit count in the case of
        :ref:`ZZ` or :ref:`QQ`.

        With ``method='GJ'`` the matrix is converted to the associated field
        domain and Gauss-Jordan elimination with division is used. With
        ``method='FF'`` fraction-free Gauss-Jordan elimination is used in the
        current domain. With ``method='CD'`` the denominators are cleared and
        fraction-free Gauss-Jordan elimination is used in the associated ring.

        If the method used is not ``'FF'`` then the computed result might be
        in a different domain than the input matrix. By default the result is
        converted back to the domain of the input matrix. To keep the result
        in whichever domain it is computed in set ``keep_domain=False``.

        By default the sparse implementations of the different algorithms are
        used. To use the dense implementations instead use e.g.
        ``method='GJ_dense'``. With ``method='auto'`` only the sparse
        implementations are used.

    keep_domain : bool, optional
        If ``True`` then the domain of the matrix is preserved. If ``False``
        then the domain might be preserved or it might be changed to an
        associated field or ring (if that is potentially faster).

    See Also
    ========

    sympy.polys.matrices.domainmatrix.DomainMatrix.rref_den
        The ``DomainMatrix`` method that calls this function.
    """
    is_dense = M.rep.fmt == 'dense'

    if method == 'auto':
        method, use_dense = _dm_rref_den_choose_method(M)
    else:
        use_dense = method.endswith('_dense')
        if use_dense:
            method = method[:-len('_dense')]

    if use_dense:
        if not is_dense:
            M = M.to_dense()

        result = _dm_rref_den(M, method, keep_domain, use_dense)
        M_rref, den, pivots = result

        if not is_dense:
            M_rref = M_rref.to_sparse()

    else:
        if is_dense:
            M = M.to_sparse()

        result = _dm_rref_den(M, method, keep_domain, use_dense)
        M_rref, den, pivots = result

        if is_dense:
            M_rref = M_rref.to_dense()

    return M_rref, den, pivots


# These are the four basic implementations that we want to choose between:


def _dm_rref_GJ_sparse(M):
    """Compute RREF using sparse Gauss-Jordan elimination with division."""
    M_rref_d, pivots, _ = sdm_irref(M.rep)
    M_rref_sdm = SDM(M_rref_d, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_sdm), pivots


def _dm_rref_GJ_dense(M):
    """Compute RREF using dense Gauss-Jordan elimination with division."""
    ddm = M.rep.copy()
    pivots = ddm_irref(ddm)
    M_rref_ddm = DDM(ddm, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_ddm), pivots


def _dm_rref_den_FF_sparse(M):
    """Compute RREF using sparse fraction-free Gauss-Jordan elimination."""
    M_rref_d, den, pivots = sdm_rref_den(M.rep, M.domain)
    M_rref_sdm = SDM(M_rref_d, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_sdm), den, pivots


def _dm_rref_den_FF_dense(M):
    """Compute RREF using sparse fraction-free Gauss-Jordan elimination."""
    ddm = M.rep.copy()
    den, pivots = ddm_irref_den(ddm, M.domain)
    M_rref_ddm = DDM(ddm, M.shape, M.domain)
    pivots = tuple(pivots)
    return M.from_rep(M_rref_ddm), den, pivots


def _dm_rref_GJ(M, use_dense):
    """Compute RREF using Gauss-Jordan with division."""
    if use_dense:
        return _dm_rref_GJ_dense(M)
    else:
        return _dm_rref_GJ_sparse(M)


def _dm_rref_den_FF(M, use_dense):
    """Compute RREF using fraction-free Gauss-Jordan."""
    if use_dense:
        return _dm_rref_den_FF_dense(M)
    else:
        return _dm_rref_den_FF_sparse(M)


def _dm_rref(M, method, use_dense):
    """Compute RREF using the given method."""

    if method == 'GJ':
        # Use Gauss-Jordan with division over the associated field.
        M_rref_f, pivots = _dm_rref_GJ(_to_field(M), use_dense)
    elif method == 'FF':
        # Use fraction-free GJ over the current domain.
        M_rref_r, den, pivots = _dm_rref_den_FF(M, use_dense)
        M_rref_f = M_rref_r.to_field() / den
    elif method == 'CD':
        # Clear denominators and use fraction-free GJ in the associated ring.
        _, Mr = M.clear_denoms(convert=True)
        M_rref_r, den, pivots = _dm_rref_den_FF(Mr, use_dense)
        M_rref_f = M_rref_r.to_field() / den
    else:
        raise ValueError(f"Unknown method for rref: {method}")

    return M_rref_f, pivots


def _dm_rref_den(M, method, keep_domain, use_dense):
    """Compute RREF with denominator using the given method."""
    K = M.domain

    if method == 'GJ':
        # Use Gauss-Jordan with division over the associated field.
        M_rref, pivots = _dm_rref_GJ(_to_field(M), use_dense)

        if keep_domain and M_rref.domain != K:
            _, M_rref = M_rref.clear_denoms(convert=True)

        if pivots:
            den = M_rref[0, pivots[0]].element
        else:
            den = K.one

    elif method == 'FF':
        # Use fraction-free GJ over the current domain.
        M_rref, den, pivots = _dm_rref_den_FF(M, use_dense)

    elif method == 'CD':
        # Clear denominators and use fraction-free GJ in the associated ring.
        _, Mr = M.clear_denoms(convert=True)
        M_rref, den, pivots = _dm_rref_den_FF(Mr, use_dense)

        if keep_domain and M_rref.domain != K:
            M_rref = M_rref.to_field() / den
            den = K.one

    else:
        raise ValueError(f"Unknown method for rref: {method}")

    return M_rref, den, pivots


def _dm_rref_choose_method(M):
    """Choose the fastest method for computing RREF for M."""

    # The sparse implementations are always faster
    use_dense = False

    K = M.domain

    if K.is_ZZ:
        method = _dm_rref_choose_method_ZZ(M)
    elif K.is_QQ:
        method = _dm_rref_choose_method_QQ(M)
    else:
        # This is definitely suboptimal. More work is needed to determine the
        # best method for computing RREF over domains that are not QQ.
        method = 'GJ'

    return method, use_dense


def _dm_rref_den_choose_method(M):
    """Choose the fastest method for computing RREF for M with denominator."""

    # The sparse implementations are always faster
    use_dense = False

    K = M.domain

    if K.is_ZZ:
        method = _dm_rref_den_choose_method_ZZ(M)
    elif K.is_QQ:
        method = _dm_rref_den_choose_method_QQ(M)
    else:
        # More work is needed to determine the best method for computing RREF
        # over domains that are not ZZ or QQ.
        method = 'FF'

    return method, use_dense


def _dm_rref_choose_method_ZZ(M):
    """Choose the fastest method for computing RREF over ZZ."""
    # Fastest method for rref is the same as for rref_den.
    return _dm_rref_den_choose_method_ZZ(M)


def _dm_rref_den_choose_method_QQ(M):
    """Choose the fastest method for computing RREF with denominator over QQ."""
    # Fastest method for rref_den is the same as for rref.
    return _dm_rref_choose_method_QQ(M)


def _dm_rref_choose_method_QQ(Mq):
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
        return 'GJ'

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
            return 'GJ'

    # If we get here then the matrix is dense and the lcm of the denominators
    # is not too large compared to the numerators. For particularly small
    # denominators it is fastest just to clear them and use fraction-free
    # Gauss-Jordan over ZZ. With very small denominators this is a little
    # faster than using rref_den over QQ but there is an intermediate regime
    # where rref_den over QQ is significantly faster. The small denominator
    # case is probably very common because small fractions like 1/2 or 1/3 are
    # often seen in user inputs.

    if denom_lcm.bit_length() < 50:
        return 'CD'
    else:
        return 'FF'


def _dm_rref_den_choose_method_ZZ(M):
    """Choose the fastest method for computing RREF with denominator over ZZ."""
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
    density, nrows_nz, ncols = _dm_row_density(M)

    # For small matrices use QQ if more than half the entries are zero.
    if nrows_nz < 10:
        if density < ncols/2:
            return 'GJ'
        else:
            return 'FF'

    # These are just shortcuts for the formula below.
    if density < 5:
        return 'GJ'
    elif density > 5 + PARAM/nrows_nz:
        return 'FF'  # pragma: no cover

    # Maximum bitsize of any entry.
    elements = _dm_elements(M)
    bits = max([e.bit_length() for e in elements], default=1)

    # Wideness parameter. This is 1 for square or tall matrices but >1 for wide
    # matrices.
    wideness = max(1, 2/3*ncols/nrows_nz)

    max_density = (5 + PARAM/(nrows_nz*bits**2)) * wideness

    if density < max_density:
        return 'GJ'
    else:
        return 'FF'


def _dm_row_density(M):
    """Density measure for sparse matrices.

    Defines the "density", ``d`` as the average number of non-zero entries per
    row except ignoring rows that are fully zero. RREF can ignore fully zero
    rows so they are excluded. By definition ``d >= 1`` except that we define
    ``d = 0`` for the zero matrix.

    Returns ``(density, nrows_nz, ncols)`` where ``nrows_nz`` counts the number
    of nonzero rows and ``ncols`` is the number of columns.
    """
    # Uses the SDM dict-of-dicts representation.
    ncols = M.shape[1]
    rows_nz = M.rep.to_sdm().values()
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


def _to_field(M):
    """Convert a DomainMatrix to a field if possible."""
    K = M.domain
    if K.has_assoc_Field:
        return M.to_field()
    else:
        return M
