from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.combinatorics.permutations import Permutation
from sympy.external.gmpy import gcd, lcm, sqrt as isqrt
from sympy.matrices.dense import MutableDenseMatrix
from sympy.ntheory import sqrt_mod, nextprime, primitive_root, primefactors
from sympy.polys.domains import ZZ, QQ, FiniteField
from sympy.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.polyclasses import ANP

if TYPE_CHECKING:
    from typing import Protocol, Iterable, Iterator, Generator, Sequence
    from sympy.combinatorics.perm_groups import PermutationGroup
    from sympy.external.gmpy import MPZ
    from sympy.polys.domains.domain import Domain

    class CC(Protocol):
        """Protocol for a conjugacy class."""

        def __len__(self) -> int: ...
        def __iter__(self) -> Iterator[Permutation]: ...
        def __contains__(self, i: Permutation, /) -> bool: ...


class CharacterTable(MutableDenseMatrix):
    r"""
    A class for creating the character table of a finite group.

    Explanation
    ============

    The character table is a complex matrix with rows corresponding to characters
    and columns corresponding to conjugacy classes. Each entry is the value of the
    character on the conjugacy class, which is the trace of the representation matrix.

    The trivial character is always the first row, and the identity class is the first column.

    Examples
    ========

    >>> from sympy.combinatorics import CharacterTable, SymmetricGroup, AlternatingGroup
    >>> CharacterTable.from_perm_group(SymmetricGroup(4)) # doctest: +SKIP
    Matrix([
     [1,  1,  1,  1,  1],
     [1, -1,  1, -1,  1],
     [2,  0,  2,  0, -1],
     [3, -1, -1,  1,  0],
     [3,  1, -1, -1,  0]])

    The character table builds upon DomainMatrix, from which users can access the
    the values of the characters as domain elements. The domain of a character table is
    either ZZ or a cyclotomic field.

    >>> M = CharacterTable.from_perm_group(AlternatingGroup(4))
    >>> M # doctest: +SKIP
    Matrix([
     [1,          1,          1,  1],
     [1, -1 - zeta3,      zeta3,  1],
     [1,      zeta3, -1 - zeta3,  1],
     [3,          0,          0, -1]])
    >>> M._rep.domain
    QQ<zeta3>
    >>> M.zeta_order
    3

    The order of the columns matches the order of the conjugacy classes,
    which can be accessed with the method `conjugacy_class_reps`.

    >>> M.conjugacy_class_reps() # doctest: +SKIP
    [Permutation(3), Permutation(0, 2, 3), Permutation(0, 1, 3), Permutation(0, 2)(1, 3)]

    References
    ==========
    .. [1] Holt, D., Eick, B., O'Brien, E.
           "Handbook of Computational Group Theory"

    .. [2] https://en.wikipedia.org/wiki/Character_table

    .. [3] https://docs.gap-system.org/doc/ref/manual.pdf

    """

    _conjugacy_class_reps: list[Permutation]

    def conjugacy_class_reps(self) -> list[Permutation]:
        return self._conjugacy_class_reps

    @property
    def zeta_order(self) -> int:
        if self._rep.domain.is_CyclotomicField:
            return self._rep.domain.zeta_order
        return 1

    @classmethod
    def from_perm_group(cls, G: PermutationGroup) -> CharacterTable:
        return dixon_character_table(G.conjugacy_classes())


def _compute_cmmatrices(
    cc: Sequence[CC], dom: Domain
) -> Generator[DomainMatrix, None, None]:
    n = len(cc)
    rmul = Permutation.rmul_with_af
    for ind in range(n):
        m = [[-1] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if m[i][j] != -1:
                    continue
                m[i] = [0] * n
                for g in cc[ind]:
                    for t in range(n):
                        out = next(iter(cc[t]))
                        # out = g^{-1} * h
                        h = rmul(out, g)
                        if h in cc[i]:
                            m[i][t] += 1
        yield DomainMatrix.from_list(m, dom)


def dixon_prime(order: int | MPZ, exponent: int | MPZ) -> int | MPZ:
    """Find a prime p so that `p > 2*sqrt(order)` and `p%exponent == 1`."""
    if order == 1:
        # trivial group => exponent == 1
        return 3
    p = 2 * isqrt(order)
    while True:
        p = nextprime(p + 1)
        if p % exponent == 1:
            break
    return p


def _simultaneous_diagonalize(
    mats: Iterable[DomainMatrix], check_dim: bool = True
) -> DomainMatrix:
    """
    Compute the simultaneous diagonalization of a list of DomainMatrices.
    Returns Z such that `Z * N * Z.inv()` is diagonal for all N in mats
    and Z is on the same domain as mats (assuming Z exists).
    """
    it = iter(mats)
    N0 = next(it)
    spaces = [S for _, __, S in N0.ground_eigenvects()]
    n = N0.shape[0]
    for N in it:
        if len(spaces) == n:
            break

        new_spaces = []
        for S in spaces:
            if S.shape[0] <= 1:
                new_spaces.append(S)
                continue
            _, pivots = S.rref()
            N0 = N.extract(range(S.shape[1]), pivots)
            for _, __, sub_S in (S * N0).transpose().ground_eigenvects():
                # sub_S: (sub_dim x dim), S: (dim x n) -> (sub_dim x n)
                sub_S = sub_S.transpose().rref()[0]
                new_spaces.append(sub_S * S)
        spaces = new_spaces

    if check_dim and len(spaces) != n:
        # not expected to happen here
        raise ValueError("Failed to compute the common eigenspace decomposition.")

    return DomainMatrix.vstack(*spaces)


def _get_invmap(cc: Sequence[CC]) -> list[int]:
    """Compute `inv_map[i] = j` such that `class[i]**-1 == class[j]`."""
    inv_map = [-1] * len(cc)
    reps = [next(iter(c)) for c in cc]
    inv_reps = [~g for g in reps]
    for i, ig in enumerate(inv_reps):
        if inv_map[i] != -1:
            continue
        for j, c in enumerate(cc):
            if ig in c:
                inv_map[i] = j
                inv_map[j] = i
                break
    return inv_map


def _get_powermap(cc: Sequence[CC], exponent: int | MPZ) -> list[list[int]]:
    """Compute `pm[i][pow] = k` such that `class[i]**pow == k`."""
    n = len(cc)
    pm = [[-1] * exponent for _ in range(n)]
    reps = [next(iter(c)) for c in cc]
    rmul = Permutation.rmul_with_af
    identity = rmul(reps[0], ~reps[0])
    for i in range(n):
        g = reps[i]
        gk = identity
        for k in range(exponent):
            for t in range(n):
                if gk in cc[t]:
                    pm[i][k] = t
                    break
            gk = rmul(gk, g)
    return pm


def _normalize_fp(cc: Sequence[CC], esd: DomainMatrix) -> list[list]:
    Fp = esd.domain
    p = Fp.mod  # type: ignore
    n = len(cc)
    rows = esd.to_list()
    cc_sizes = [Fp(len(c)) for c in cc]
    G_order = sum(len(cc[t]) for t in range(len(cc)))

    inv_map = _get_invmap(cc)
    normalized_rows = []
    for row in rows:
        # 1. normalize so the first class is 1
        scale = 1 / Fp(row[0])
        row = [x * scale for x in row]

        # 2. chi(1)^2 * sum |Ci| * row[i] * row[inv_i] = |G|
        dot = sum(cc_sizes[k] * row[k] * row[inv_map[k]] for k in range(n))

        # chi1_sq = |G| / dot
        chi1_sq = Fp(G_order) / dot

        root = sqrt_mod(int(chi1_sq), p)
        if root is None:
            # not expected to happen
            raise ValueError(f"Failed to compute sqrt({chi1_sq}) on {Fp}")
        chi1 = Fp(root)

        normalized_rows.append([x * chi1 for x in row])

    return normalized_rows


def _get_global_conductor(
    rows: list[list], pm: list[list[int]], m: int | MPZ
) -> int | MPZ:
    """
    Find the smallest natural number k such that all values of the
    character table can be embedded in the k-th cyclotomic field.
    """
    n = len(rows)
    gal_m = [a for a in range(m) if gcd(a, m) == 1]
    p_factors = primefactors(m)[::-1]

    for p in p_factors:
        while m % p == 0:
            # test invariance under m//p
            m = m // p
            r = 1 % m

            fixed = True
            for a in gal_m:
                if a % m != r:
                    continue
                for i in range(n):
                    row = rows[i]
                    for j in range(n):
                        if row[pm[j][a]] != row[j]:
                            fixed = False
                            break
                    if not fixed:
                        break
                if not fixed:
                    break

            if not fixed:
                m = m * p  # restore
                break
    return m


def _lift_to_minimal_field(normalized_rows, pm, k, e, Fp):
    n = len(normalized_rows)
    p = Fp.mod
    half_p = p // 2

    if k == 1:
        # integer character table
        int_rows = []
        for i in range(n):
            # abs(every entry) <= sqrt(|G|) <= p//2
            row = [ZZ(int(v)) for v in normalized_rows[i]]
            row = [v if v <= half_p else v - p for v in row]
            int_rows.append(row)
        _sort_characters(int_rows, ZZ)
        return DomainMatrix(int_rows, (n, n), ZZ)

    dom = QQ.cyclotomic_field(k)

    x = Fp(primitive_root(p)) ** ((p - 1) // k)

    gal = [a for a in range(k) if gcd(a, k) == 1]
    phi = len(gal)

    # V[a, i] = (x**a)**i
    V = []
    for a in gal:
        xa = x**a
        row = [xa**i for i in range(phi)]
        V.append(row)

    V_mat = DomainMatrix(V, (phi, phi), Fp)
    V_inv = V_mat.inv()

    # for a in (Z/k)*, find A in (Z/exp)* that A = a (mod k)
    gal_lifted = []
    for a in gal:
        A = a
        while gcd(A, e) != 1:
            A += k
        gal_lifted.append(A)

    dM = []
    for i in range(n):
        char_row_fp = normalized_rows[i]
        row_anps = []

        for j in range(n):
            # [chi(g^A1), chi(g^A2), ...] mod p
            b_vec_data = []
            for A in gal_lifted:
                class_of_gA = pm[j][A % e]
                b_vec_data.append([char_row_fp[class_of_gA]])

            b_vec = DomainMatrix(b_vec_data, (phi, 1), Fp)

            # c = V^-1 * b
            c_vec = (V_inv * b_vec).to_list_flat()

            c_ints = [int(v) for v in c_vec]
            c_ints = [v if v <= half_p else v - p for v in c_ints]

            row_anps.append(ANP(c_ints[::-1], dom.mod, QQ))

        dM.append(row_anps)

    _sort_characters(dM, dom)
    return DomainMatrix(dM, (n, n), dom)


def _sort_characters(rows: list[list], dom: Domain):
    """
    Sort the character table and move the trivial
    character to the first row. Done in-place.
    """
    one = dom.one
    rows.sort()
    for i in range(len(rows)):
        if all(v == one for v in rows[i]):
            rows[0], rows[i] = rows[i], rows[0]
            break


def dixon_character_table(conjugacy_classes: Sequence[CC]) -> CharacterTable:
    """
    Compute the character table of a finite group from its conjugacy classes
    using Dixon's algorithm.

    Parameters
    ==========
    conjugacy_classes : Sequence[CC]
        The conjugacy classes of the group.

    References
    ==========
    .. [1] Dixon, J.
             "High Speed Computation of Group Characters"

    .. [2] Schneider, J.
             "Dixon's character table algorithm revisited""

    .. [3] Holt, D., Eick, B., O'Brien, E.
             "Handbook of Computational Group Theory"
    """
    cc = conjugacy_classes
    order = sum(len(cc[t]) for t in range(len(cc)))
    exponent = lcm(*(int(next(iter(c)).order()) for c in cc))
    p = dixon_prime(order, exponent)
    Fp = FiniteField(p)

    mats = _compute_cmmatrices(cc, Fp)
    esd = _simultaneous_diagonalize(mats)

    normalized = _normalize_fp(cc, esd)
    pm = _get_powermap(cc, exponent)

    conductor = _get_global_conductor(normalized, pm, exponent)
    dm = _lift_to_minimal_field(normalized, pm, conductor, exponent, Fp)

    tbl = CharacterTable._fromrep(dm)
    tbl._conjugacy_class_reps = [next(iter(c)) for c in cc]
    return tbl
