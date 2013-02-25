"Implementation of matrix FGLM Groebner basis conversion algorithm. """

from sympy.polys.distributedpolys import sdp_strip, sdp_sort, sdp_sub, sdp_rem, sdp_monic, sdp_LM
from sympy.polys.monomialtools import monomial_mul, monomial_div

def matrix_fglm(F, u, O_from, O_to, K):
    """
    Converts the reduced Groebner basis ``F`` of a zero-dimensional
    ideal w.r.t. ``O_from`` to a reduced Groebner basis
    w.r.t. ``O_to``.

    References
    ==========

    J.C. Faugere, P. Gianni, D. Lazard, T. Mora (1994). Efficient
    Computation of Zero-dimensional Groebner Bases by Change of
    Ordering
    """
    old_basis = _basis(F, u, O_from, K)
    M = _representing_matrices(old_basis, F, u, O_from, K)

    # V contains the normalforms (wrt O_from) of S
    S = [(0,) * (u + 1)]
    V = [[K.one] + [K.zero] * (len(old_basis) - 1)]
    G = []

    L = [(i, 0) for i in xrange(u + 1)]  # (i, j) corresponds to x_i * S[j]
    L.sort(key=lambda (k, l): O_to(_incr_k(S[l], k)), reverse=True)
    t = L.pop()

    P = _identity_matrix(len(old_basis), K)

    while True:
        s = len(S)
        v = _matrix_mul(M[t[0]], V[t[1]], K)
        _lambda = _matrix_mul(P, v, K)

        if all(_lambda[i] == K.zero for i in xrange(s, len(old_basis))):
            # there is a linear combination of v by V

            lt = [(_incr_k(S[t[1]], t[0]), K.one)]
            rest = sdp_strip(sdp_sort([(S[i], _lambda[i]) for i in xrange(s)], O_to))
            g = sdp_sub(lt, rest, u, O_to, K)

            if g != []:
                G.append(g)

        else:
            # v is linearly independant from V
            P = _update(s, _lambda, P, K)
            S.append(_incr_k(S[t[1]], t[0]))
            V.append(v)

            L.extend([(i, s) for i in xrange(u + 1)])
            L = list(set(L))
            L.sort(key=lambda (k, l): O_to(_incr_k(S[l], k)), reverse=True)

        L = [(k, l) for (k, l) in L if
            all(monomial_div(_incr_k(S[l], k), sdp_LM(g, u)) is None for g in G)]

        if not L:
            G = [ sdp_monic(g, K) for g in G ]
            return sorted(G, key=lambda g: O_to(sdp_LM(g, u)), reverse=True)

        t = L.pop()


def _incr_k(m, k):
    return tuple(list(m[:k]) + [m[k] + 1] + list(m[k + 1:]))


def _identity_matrix(n, K):
    M = [[K.zero] * n for _ in xrange(n)]

    for i in xrange(n):
        M[i][i] = K.one

    return M


def _matrix_mul(M, v, K):
    return [sum([row[i] * v[i] for i in xrange(len(v))]) for row in M]


def _update(s, _lambda, P, K):
    """
    Update ``P`` such that for the updated `P'` `P' v = e_{s}`.
    """
    k = min([j for j in xrange(s, len(_lambda)) if _lambda[j] != 0])

    for r in xrange(len(_lambda)):
        if r != k:
            P[r] = [P[r][j] - (
                P[k][j] * _lambda[r]) / _lambda[k] for j in xrange(len(P[r]))]

    P[k] = [P[k][j] / _lambda[k] for j in xrange(len(P[k]))]

    P[k], P[s] = P[s], P[k]

    return P


def _representing_matrices(basis, G, u, O, K):
    """
    Compute the matrices corresponding to the linear maps `m \mapsto
    x_i m` for all variables `x_i`.
    """
    def var(i):
        return tuple([0] * i + [1] + [0] * (u - i))

    def representing_matrix(m):
        M = [[K.zero] * len(basis) for _ in xrange(len(basis))]

        for i, v in enumerate(basis):
            r = sdp_rem([(monomial_mul(m, v), K.one)], G, u, O, K)

            for term in r:
                j = basis.index(term[0])
                M[j][i] = term[1]

        return M

    return [representing_matrix(var(i)) for i in xrange(u + 1)]


def _basis(G, u, O, K):
    """
    Computes a list of monomials which are not divisible by the leading
    monomials wrt to ``O`` of ``G``. These monomials are a basis of
    `K[X_1, \ldots, X_n]/(G)`.
    """
    leading_monomials = [sdp_LM(g, u) for g in G]
    candidates = [(0,) * (u + 1)]
    basis = []

    while candidates:
        t = candidates.pop()
        basis.append(t)

        new_candidates = [_incr_k(t, k) for k in xrange(u + 1)
            if all(monomial_div(_incr_k(t, k), lmg) is None
            for lmg in leading_monomials)]
        candidates.extend(new_candidates)
        candidates.sort(key=lambda m: O(m), reverse=True)

    basis = list(set(basis))

    return sorted(basis, key=lambda m: O(m))
