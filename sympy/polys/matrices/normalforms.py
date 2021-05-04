'''Functions returning normal forms of matrices'''

from .domainmatrix import DomainMatrix


def smith_normal_form(m):
    '''
    Return the Smith Normal Form of a matrix `m` over the ring `domain`.
    This will only work if the ring is a principal ideal domain.

    Examples
    ========

    >>> from sympy import ZZ
    >>> from sympy.polys.matrices import DomainMatrix
    >>> from sympy.polys.matrices.normalforms import smith_normal_form
    >>> m = DomainMatrix([[ZZ(12), ZZ(6), ZZ(4)],
    ...                   [ZZ(3), ZZ(9), ZZ(6)],
    ...                   [ZZ(2), ZZ(16), ZZ(14)]], (3, 3), ZZ)
    >>> print(smith_normal_form(m).to_Matrix())
    Matrix([[1, 0, 0], [0, 10, 0], [0, 0, -30]])

    '''
    invs = invariant_factors(m)
    smf = DomainMatrix.diag(invs, m.domain, m.shape)
    return smf


def invariant_factors(m):
    '''
    Return the tuple of abelian invariants for a matrix `m`
    (as in the Smith-Normal form)

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Smith_normal_form#Algorithm
    [2] http://sierra.nmsu.edu/morandi/notes/SmithNormalForm.pdf

    '''
    domain = m.domain
    if not domain.is_PID:
        msg = "The matrix entries must be over a principal ideal domain"
        raise ValueError(msg)

    if 0 in m.shape:
        return ()

    rows, cols = shape = m.shape
    m = list(m.to_dense().rep)

    def add_rows(m, i, j, a, b, c, d):
        # replace m[i, :] by a*m[i, :] + b*m[j, :]
        # and m[j, :] by c*m[i, :] + d*m[j, :]
        for k in range(cols):
            e = m[i][k]
            m[i][k] = a*e + b*m[j][k]
            m[j][k] = c*e + d*m[j][k]

    def add_columns(m, i, j, a, b, c, d):
        # replace m[:, i] by a*m[:, i] + b*m[:, j]
        # and m[:, j] by c*m[:, i] + d*m[:, j]
        for k in range(rows):
            e = m[k][i]
            m[k][i] = a*e + b*m[k][j]
            m[k][j] = c*e + d*m[k][j]

    def clear_column(m):
        # make m[1:, 0] zero by row and column operations
        if m[0][0] == 0:
            return m  # pragma: nocover
        pivot = m[0][0]
        for j in range(1, rows):
            if m[j][0] == 0:
                continue
            d, r = domain.div(m[j][0], pivot)
            if r == 0:
                add_rows(m, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[j][0])
                d_0 = domain.div(m[j][0], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_rows(m, 0, j, a, b, d_0, -d_j)
                pivot = g
        return m

    def clear_row(m):
        # make m[0, 1:] zero by row and column operations
        if m[0][0] == 0:
            return m  # pragma: nocover
        pivot = m[0][0]
        for j in range(1, cols):
            if m[0][j] == 0:
                continue
            d, r = domain.div(m[0][j], pivot)
            if r == 0:
                add_columns(m, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[0][j])
                d_0 = domain.div(m[0][j], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_columns(m, 0, j, a, b, d_0, -d_j)
                pivot = g
        return m

    # permute the rows and columns until m[0,0] is non-zero if possible
    ind = [i for i in range(rows) if m[i][0] != 0]
    if ind and ind[0] != 0:
        m[0], m[ind[0]] = m[ind[0]], m[0]
    else:
        ind = [j for j in range(cols) if m[0][j] != 0]
        if ind and ind[0] != 0:
            for row in m:
                row[0], row[ind[0]] = row[ind[0]], row[0]

    # make the first row and column except m[0,0] zero
    while (any([m[0][i] != 0 for i in range(1,cols)]) or
           any([m[i][0] != 0 for i in range(1,rows)])):
        m = clear_column(m)
        m = clear_row(m)

    if 1 in shape:
        invs = ()
    else:
        lower_right = DomainMatrix([r[1:] for r in m[1:]], (rows-1, cols-1), domain)
        invs = invariant_factors(lower_right)

    if m[0][0]:
        result = [m[0][0]]
        result.extend(invs)
        # in case m[0] doesn't divide the invariants of the rest of the matrix
        for i in range(len(result)-1):
            if result[i] and domain.div(result[i+1], result[i])[1] != 0:
                g = domain.gcd(result[i+1], result[i])
                result[i+1] = domain.div(result[i], g)[0]*result[i+1]
                result[i] = g
            else:
                break
    else:
        result = invs + (m[0][0],)
    return tuple(result)
