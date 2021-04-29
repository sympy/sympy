'''Functions returning normal forms of matrices'''

def smith_normal_form(m, domain=None, full=False):
    '''
    Return the Smith Normal Form of a matrix ``m`` over the ring ``domain``.
    This will only work if the ring is a principal ideal domain.
    Output is a matrix of the same type as ``m``.

    Examples
    ========

    >>> from sympy.polys.polymatrix import PolyMatrix as Matrix
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.matrices.normalforms import smith_normal_form
    >>> m = Matrix([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
    >>> setattr(m, "ring", ZZ)
    >>> print(smith_normal_form(m))
    Matrix([[mpz(1), 0, 0], [0, mpz(10), 0], [0, 0, mpz(-30)]])

    Output when ``full=True`` is a tuple with the element the same as i
    ``full=False`` and next two elements will be the invertible matrices
    ``s, t`` such that the product ``s, m, t`` is the Smith Normal Form.

    >>> smf, s, t = smith_normal_form(m, full=True)
    >>> print(s)
    Matrix([[0, 1, -1], [1, -4, 0], [0, 2, -3]])
    >>> print(t)
    Matrix([[1, 1, 10], [0, -1, -2], [0, 1, 3]])

    '''
    if full:
        invs, s, t = invariant_factors(m, domain=domain, full=True)
    else:
        invs = invariant_factors(m, domain=domain, full=False)
    smf = m.diag(*invs)
    n = len(invs)
    if m.rows > n:
        smf = smf.row_insert(m.rows, m.zeros(m.rows-n, m.cols))
    elif m.cols > n:
        smf = smf.col_insert(m.cols, m.zeros(m.rows, m.cols-n))
    if full:
        return smf, s, t
    return smf

def invariant_factors(m, domain=None, full=False):
    '''
    Return the tuple of abelian invariants for a matrix ``m``
    (as in the Smith-Normal form)

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Smith_normal_form#Algorithm
    [2] http://sierra.nmsu.edu/morandi/notes/SmithNormalForm.pdf

    '''
    if not domain:
        if not (hasattr(m, "ring") and m.ring.is_PID):
            raise ValueError(
                "The matrix entries must be over a principal ideal domain")
        else:
            domain = m.ring

    if len(m) == 0:
        if full:
            return (), m, m
        return ()

    m = m[:, :]
    if full:
        s = m.eye(m.rows)
        t = m.eye(m.cols)
    else:
        s = None
        t = None

    def add_rows(m, s, t, i, j, a, b, c, d):
        # replace m[i, :] by a*m[i, :] + b*m[j, :]
        # and m[j, :] by c*m[i, :] + d*m[j, :]
        mats = (m,)
        if full:
            mats = (m, s)
        for mat in mats:
            for k in range(mat.cols):
                e = mat[i, k]
                mat[i, k] = a*e + b*mat[j, k]
                mat[j, k] = c*e + d*mat[j, k]

    def add_columns(m, s, t, i, j, a, b, c, d):
        # replace m[:, i] by a*m[:, i] + b*m[:, j]
        # and m[:, j] by c*m[:, i] + d*m[:, j]
        mats = (m,)
        if full:
            mats = (m, t)
        for mat in mats:
            for k in range(mat.rows):
                e = mat[k, i]
                mat[k, i] = a*e + b*mat[k, j]
                mat[k, j] = c*e + d*mat[k, j]

    def clear_column(m, s, t):
        # make m[1:, 0] zero by row and column operations
        if m[0,0] == 0:
            return
        pivot = m[0, 0]
        for j in range(1, m.rows):
            if m[j, 0] == 0:
                continue
            d, r = domain.div(m[j,0], pivot)
            if r == 0:
                add_rows(m, s, t, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[j,0])
                d_0 = domain.div(m[j, 0], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_rows(m, s, t, 0, j, a, b, d_0, -d_j)
                pivot = g

    def clear_row(m, s, t):
        # make m[0, 1:] zero by row and column operations
        if m[0] == 0:
            return
        pivot = m[0, 0]
        for j in range(1, m.cols):
            if m[0, j] == 0:
                continue
            d, r = domain.div(m[0, j], pivot)
            if r == 0:
                add_columns(m, s, t, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[0, j])
                d_0 = domain.div(m[0, j], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_columns(m, s, t, 0, j, a, b, d_0, -d_j)
                pivot = g

    # permute the rows and columns until m[0,0] is non-zero if possible
    ind = [i for i in range(m.rows) if m[i,0] != 0]
    if ind and ind[0] != 0:
        m = m.permute_rows([[0, ind[0]]])
        if full:
            s = s.permute_rows([[0, ind[0]]])
    else:
        ind = [j for j in range(m.cols) if m[0,j] != 0]
        if ind and ind[0] != 0:
            m = m.permute_cols([[0, ind[0]]])
            if full:
                t = t.permute_cols([[0, ind[0]]])

    # make the first row and column except m[0,0] zero
    while (any([m[0,i] != 0 for i in range(1,m.cols)]) or
           any([m[i,0] != 0 for i in range(1,m.rows)])):
        clear_column(m, s, t)
        clear_row(m, s, t)

    if 1 in m.shape:
        invs = ()
    else:
        if full:
            invs, s_small, t_small = invariant_factors(m[1:,1:], domain=domain,
                                        full=True)
            s2 = m.eye(m.rows)
            t2 = m.eye(m.cols)
            s2[1:,1:] = s_small
            t2[1:,1:] = t_small
            s = s2 * s
            t = t * t2
        else:
            invs = invariant_factors(m[1:,1:], domain=domain, full=False)

    if m[0,0]:
        result = [m[0,0]]
        result.extend(invs)
        # in case m[0] doesn't divide the invariants of the rest of the matrix
        for i in range(len(result)-1):
            if result[i] and domain.div(result[i+1], result[i])[1] != 0:
                g = domain.gcd(result[i+1], result[i])
                p = domain.div(result[i], g)[0]
                if full:
                    s[i+1, :] *= p
                    s[i, :] *= g/result[i]
                result[i+1] *= p
                result[i] = g
            else:
                break
    else:
        if full:
            if s.rows > 1:
                permute_list = [[i, i+1] for i in range(s.rows-1)]
                s = s.permute_rows(permute_list)
            if t.cols > 1:
                permute_list = [[i, i+1] for i in range(t.cols-1)]
                t = t.permute_cols(permute_list)
        result = invs + (m[0,0],)
    if full:
        return tuple(result), s, t
    return tuple(result)
