from sympy.matrices import *

'''
Functions returning normal forms of matrices

'''
def smith_normal_form(m, domain = None):
    invs = smith_normal_invariants(m)
    smf = diag(*invs)
    n = len(invs)
    if m.rows > n:
        smf = smf.row_insert(m.rows, zeros(m.rows-n, m.cols))
    elif m.cols > n:
        smf = smf.col_insert(m.cols, zeros(m.rows, m.cols-n))
    return smf

def smith_normal_invariants(m, domain = None):
    '''
    Return the tuple of abelian invariants for a matrix `m`
    (as in the Smith-Normal form)

    '''
    if not domain and not (hasattr(m, "ring") and m.ring.is_PID):
        raise ValueError(
            "The matrix entries must be over a principal ideal domain")
    elif not domain:
        domain = m.ring

    m = m[:, :]

    def add_rows(m, i, j, a, b, c, d):
        # replace m[i, :] by a*m[i, :] + b*m[j, :]
        # and m[j, :] by c*m[i, :] + d*m[j, :]
        for k in range(m.cols):
            e = m[i, k]
            m[i, k] = a*e + b*m[j, k]
            m[j, k] = c*e + d*m[j, k]

    def add_columns(m, i, j, a, b, c, d):
        # replace m[:, i] by a*m[:, i] + b*m[:, j]
        # and m[:, j] by c*m[:, i] + d*m[:, j]
        for k in range(m.rows):
            e = m[k, i]
            m[k, i] = a*e + b*m[k, j]
            m[k, j] = c*e + d*m[k, j]

    def clear_column(m):
        # make m[1:, 0] zero by row and column operations
        ind = [i for i in range(m.rows) if m[i,0] != 0]
        if ind:
            m = m.permute_rows([[0, ind[0]]])
        else:
            return m
        pivot = m[0, 0]
        for j in range(1, m.rows):
            if m[j, 0] == 0:
                continue
            d, r = domain.div(m[j,0], pivot)
            if r == 0:
                add_rows(m, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[j,0])
                d_0 = domain.div(m[j, 0], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_rows(m, 0, j, a, b, d_0, -d_j)
                pivot = g
        return m

    def clear_row(m):
        # make m[0, 1:] zero by row and column operations
        ind = [j for j in range(m.cols) if m[0,j] != 0]
        if ind:
            m = m.permute_cols([[0, ind[0]]])
        else:
            return m
        pivot = m[0, 0]
        for j in range(1, m.cols):
            if m[0, j] == 0:
                continue
            d, r = domain.div(m[0, j], pivot)
            if r == 0:
                add_columns(m, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[0, j])
                d_0 = domain.div(m[0, j], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_columns(m, 0, j, a, b, d_0, -d_j)
                pivot = g
        return m

    m = clear_column(m)
    m = clear_row(m)
    if 1 in m.shape:
        return [m[0]]

    while not (m[0,1:].is_zero) or not (m[1:,0].is_zero):
        m = clear_column(m)
        m = clear_row(m)

    invs = smith_normal_invariants(m[1:,1:], domain=domain)
    zeros = []
    result = []
    for r in invs:
        if r == 0:
            zeros.append(r)
        else:
            result.append(r)
    new = m[0]
    # in case m[0] doesn't divide the invariants of the rest of the matrix
    if result and new != 0 and domain.div(result[0], new)[1] != 0:
        g = domain.gcd(new, result[0])
        result[0] = domain.div(new, g)[0]*result[0]
        new = g
    result = zeros + [new] + result
    return tuple(result)
