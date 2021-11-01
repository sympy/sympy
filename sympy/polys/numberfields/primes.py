"""Prime ideals in number fields. """

from sympy.core.symbol import symbols
from sympy.polys.polytools import Poly
from sympy.polys.domains.finitefield import FF
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.domains.integerring import ZZ
from sympy.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.polyutils import IntegerPowerable
from sympy.utilities.decorator import public
from .basis import round_two, nilradical_mod_p
from .modules import ModuleEndomorphism, Order, find_min_poly
from .utilities import coeff_search, supplement_a_subspace


class PrimeIdeal(IntegerPowerable):
    """
    A prime ideal in a ring of algebraic integers.
    """

    def __init__(self, ZK, p, alpha, f, e=None):
        """
        Parameters
        ==========

        ZK: the :py:class:`Order` where this ideal lives
        p: the rational prime this ideal divides.
        alpha: :py:class:`PowerBasisElement` such that the ideal is equal
            to ``p*ZK + alpha*ZK``.
        f: the inertia degree.
        e: the ramification index, if already known. If ``None``, we will
            compute it here.

        """
        self.ZK = ZK
        self.p = p
        self.alpha = alpha
        self.f = f
        self._beta = None
        self.e = e if e is not None else self.valuation(p * ZK)

    def pretty(self, theta=None):
        theta = theta or symbols('x')
        p, alpha, e, f = self.p, self.alpha, self.e, self.f
        alpha_rep = str(alpha.numerator(x=theta).as_expr())
        if alpha.denom > 1:
            alpha_rep = f'({alpha_rep})/{alpha.denom}'
        return f'[ ({p}, {alpha_rep}) e={e}, f={f} ]'

    def __repr__(self):
        return self.pretty()

    def as_submodule(self):
        return self.p * self.ZK + self.alpha * self.ZK

    def copy(self):
        return PrimeIdeal(self.ZK.copy(), self.p, self.alpha.copy(), self.f, self.e)

    def __eq__(self, other):
        if isinstance(other, PrimeIdeal):
            return self.as_submodule() == other.as_submodule()
        return NotImplemented

    def __mul__(self, other):
        return self.as_submodule() * other

    __rmul__ = __mul__

    def _zeroth_power(self):
        return self.ZK.copy()

    def _first_power(self):
        return self.copy()

    def beta(self):
        """
        Write P for this prime ideal, p for the rational prime it divides.
        Then, for computing P-adic valuations it is useful to have a number
        beta in ZK such that p/P = p*ZK + beta*ZK. Essentially, this is the
        number Psi (or the "reagent") from Kummer's 1847 paper (*Ueber die
        Zerlegung...*, Crelle vol. 35).
        """
        if self._beta is None:
            self._beta = _compute_beta(self.p, [self.alpha], self.ZK)
        return self._beta

    def valuation(self, I):
        """
        Compute the P-adic valuation of integral ideal I, given as an HNF.
        """
        return prime_valuation(I, self)


def _compute_beta(p, gens, ZK):
    r"""
    Compute the "beta" for a :ref:`PrimeIdeal` $P$ (see ``beta`` method in
    that class).

    Parameters
    ==========

    p: the rational prime $P$ divides

    gens: list of :py:class:`PowerBasisElement`s being a complete set of
        generators for $P$ over *ZK*, EXCEPT that an element equivalent to
        rational *p* can and should be omitted (since it has no effect except
        to waste time).

    ZK: maximal :py:class:`Order`

    Returns
    =======

    PowerBasisElement for beta

    References
    ==========

    [1] Cohen, H. *A Course in Computational Algebraic Number Theory.*
    (See Proposition 4.8.15.)

    """
    E = ZK.endomorphism_ring()
    matrices = [E.inner_endomorphism(g).matrix(modulus=p) for g in gens]
    B = DomainMatrix.zeros((0, ZK.n), FF(p)).vstack(*matrices)
    # A nonzero element of the nullspace of B will represent a
    # lin comb over the omegas which (i) is not a multiple of p
    # (since it is nonzero over FF(p)), while (ii) is such that
    # its product with each g in gens _is_ a multiple of p (since
    # B represents multiplication by these generators). Theory
    # predicts that such an element must exist, so nullspace should
    # be non-trivial.
    x = B.nullspace()[0, :].transpose()
    beta = ZK.pb_elt_from_col(ZK.matrix * x, denom=ZK.denom)
    return beta


@public
def prime_valuation(I, P):
    r"""
    Compute the p-adic valuation for an integral ideal *I* at a given prime
    ideal *P*.

    Parameters
    ==========

    I: :py:class:`Ideal`, representing an integral ideal whose valuation is desired.

    P: :py:class:`PrimeIdeal` at which to compute the valuation.

    Returns
    =======

    int (non-negative since *I* is required to be an integral ideal) giving the
    p-adic valuation of *I* at the prime *P*.

    Examples
    ========

    >>> from sympy import QQ
    >>> from sympy.abc import theta
    >>> from sympy.polys import cyclotomic_poly
    >>> from sympy.polys.numberfields import prime_valuation
    >>> T = cyclotomic_poly(5)
    >>> K = QQ.algebraic_field((T, theta))
    >>> P = K.primes_above(5)
    >>> ZK = K.maximal_order()
    >>> print(prime_valuation(25*ZK, P[0]))
    8

    References
    ==========

    [1] Cohen, H. *A Course in Computational Algebraic Number Theory.*
    (See Algorithm 4.8.17.)

    """
    p, ZK = P.p, P.ZK
    n, W, d = ZK.n, ZK.matrix, ZK.denom

    A = W.convert_to(QQ).inv() * I.matrix * d / I.denom
    # Although A must have integer entries, given that I is an integral ideal,
    # as a DomainMatrix it will still be over QQ, so we convert back:
    A = A.convert_to(ZZ)
    D = A.det()
    if D % p != 0:
        return 0

    beta = P.beta()

    f = d ** n // W.det()
    need_complete_test = (f % p == 0)
    v = 0
    while True:
        # Entering the loop, the cols of A represent lin combs of omegas.
        # Turn them into lin combs of thetas:
        A = W * A
        # And then one column at a time...
        for j in range(n):
            c = ZK.pb_elt_from_col(A[:, j], denom=d)
            c *= beta
            # ...turn back into lin combs of omegas, after multiplying by beta:
            c = ZK.represent(c).flat()
            for i in range(n):
                A[i, j] = c[i]
        if A[n - 1, n - 1].element % p != 0:
            break
        A = A / p
        # As noted above, domain converts to QQ even when division goes evenly.
        # So must convert back, even when we don't "need_complete_test".
        if need_complete_test:
            # In this case, having a non-integer entry is actually just our
            # halting condition.
            try:
                A = A.convert_to(ZZ)
            except CoercionFailed:
                break
        else:
            # In this case theory says we should not have any non-integer entries.
            A = A.convert_to(ZZ)
        v += 1
    return v


def _two_elt_rep(gens, ZK, p, f=None, Np=None):
    r"""
    Given a set of *ZK*-generators of a prime ideal, compute a set of just two
    *ZK*-generators for the same ideal, one of which is *p* itself.

    Parameters
    ==========

    gens: list of :py:class:`PowerBasisElement`s giving generators for the
        prime ideal over *ZK*, the ring of integers of the field $K$.

    ZK: maximal :py:class:`Order` in $K$.

    p: the rational prime divided by the prime ideal.

    f: (optional): the inertia degree of the prime ideal, if known.

    Np: (optional): the norm $p^f$ of the prime ideal, if known.
        NOTE: There is no reason to supply both *f* and *Np*. Either one will
        save us from having to compute the norm *Np* ourselves. If both are known,
        *Np* is preferred since it saves one exponentiation.

    Returns
    =======

    PowerBasisElement representing a single algebraic integer alpha such that
    the prime ideal is equal to ``p*ZK + alpha*ZK``.

    References
    ==========

    [1] Cohen, H. *A Course in Computational Algebraic Number Theory.*
    (See Algorithm 4.7.10.)

    """
    pb = ZK.parent
    T = ZK.T
    # Detect the special cases in which either (a) all generators are multiples
    # of p, or (b) there are no generators (so `all` is vacuously true):
    if all((g % p).equiv(0) for g in gens):
        return pb.zero()

    if Np is None:
        if f is not None:
            Np = p**f
        else:
            Np = abs(pb.submodule_from_gens(gens).matrix.det())

    assert isinstance(ZK, Order)
    omega = ZK.power_basis_elements()
    beta = [p*om for om in omega[1:]]  # note: we omit omega[0] == 1
    beta += gens
    search = coeff_search(len(beta), 1)
    for c in search:
        alpha = sum(ci*betai for ci, betai in zip(c, beta))
        # Note: It may be tempting to reduce alpha mod p here, to try to work
        # with smaller numbers, but must not do that, as it can result in an
        # infinite loop! E.g. try factoring 2 in Q(sqrt(-7)).
        n = alpha.norm(T) // Np
        if n % p != 0:
            # Now can reduce alpha mod p.
            return alpha % p


def _prime_decomp_easy_case(p, ZK):
    r"""
    Compute the decomposition of rational prime *p* in the ring of integers *ZK*
    (as an :py:class:`Order`), in the "easy case", i.e. the case where *p* does
    not divide the index of $\theta$ in *ZK*, $\theta$ a root of ``ZK.T``.
    """
    T = ZK.T
    T_bar = Poly(T, modulus=p)
    lc, fl = T_bar.factor_list()
    return [PrimeIdeal(ZK, p,
                       ZK.pb_elt_from_poly(Poly(t, domain=ZZ)), t.degree(), e)
            for t, e in fl]


def _prime_decomp_compute_kernel(I, p, ZK):
    """
    Parameters
    ==========

    I: a Module, representing an ideal of ``ZK/pZK``.
    p: the rational prime being factored
    ZK: maximal Order

    Returns
    =======

    Pair ``(N, G)``, where:
      N: Module representing the kernel of the map ``a |--> a**p - a`` on ``(O/pO)/I``,
         guaranteed to be a module with unity, along with
      G: Module representing a basis for the separable algebra ``A = O/I`` (see Cohen).

    """
    W = I.matrix
    n, r = W.shape
    # Want to take the Fp-basis given by the columns of I, adjoin (1, 0, ..., 0)
    # (which we know is not already in there since I is a basis for a prime ideal)
    # and then supplement this with additional columns to make an invertible n x n
    # matrix. This will then represent a full basis for ZK, whose first r columns
    # are pullbacks of the basis for I.
    if r == 0:
        B = W.eye(n, ZZ)
    else:
        B = W.hstack(W.eye(n, ZZ)[:, 0])
    if B.shape[1] < n:
        B = supplement_a_subspace(B.convert_to(FF(p))).convert_to(ZZ)

    G = ZK.submodule_from_matrix(B)
    # Must compute G's multiplication table _before_ discarding the first r
    # columns. (See Step 9 in Alg 6.2.9 in Cohen, where the betas are actually
    # needed in order to represent each product of gammas. However, once we've
    # found the representations, then we can ignore the betas.)
    G.compute_mult_tab()
    G = G.discard_before(r)

    phi = ModuleEndomorphism(G, lambda x: x**p - x)
    N = phi.kernel(modulus=p)
    assert N.starts_with_unity()
    return N, G


def _prime_decomp_maximal_ideal(I, p, ZK):
    """
    We have reached the case where we have a maximal (hence prime) ideal *I*,
    which we know because the quotient ``O/I`` is a field.

    Parameters
    ==========

    I: a Module representing an ideal of ``O/pO``.
    p: the rational prime being factored
    ZK: maximal Order

    Returns
    =======

    PrimeIdeal instance representing this prime

    """
    m, n = I.matrix.shape
    f = m - n
    G = ZK.matrix * I.matrix
    gens = [ZK.pb_elt_from_col(G[:, j], denom=ZK.denom) for j in range(G.shape[1])]
    alpha = _two_elt_rep(gens, ZK, p, f=f)
    return PrimeIdeal(ZK, p, alpha, f)


def _prime_decomp_split_ideal(I, p, N, G, ZK):
    """
    Perform the step in the prime decomposition algorithm where we have determined
    the the quotient ``ZK/I`` is _not_ a field, and we want to perform a non-trivial
    factorization of *I* by locating an idempotent element of ``ZK/I``.
    """
    assert I.parent == ZK and G.parent is ZK and N.parent is G
    # Since ZK/I is not a field, the kernel computed in the previous step contains
    # more than just the prime field Fp, and our basis N for the nullspace therefore
    # contains at least a second column (which represents an element outside Fp).
    # Let alpha be such an element:
    alpha = N(1).to_parent()
    assert alpha.module is G

    alpha_powers = []
    m = find_min_poly(alpha, FF(p), powers=alpha_powers)
    # TODO (future work):
    #  We don't actually need full factorization, so might use a faster method
    #  to just break off a single non-constant factor m1?
    lc, fl = m.factor_list()
    m1 = fl[0][0]
    m2 = m.quo(m1)
    U, V, g = m1.gcdex(m2)
    # Sanity check: theory says m is squarefree, so m1, m2 should be coprime:
    assert g == 1
    E = list(reversed(Poly(U * m1, domain=ZZ).rep.rep))
    eps1 = sum(E[i]*alpha_powers[i] for i in range(len(E)))
    eps2 = 1 - eps1
    idemps = [eps1, eps2]
    factors = []
    for eps in idemps:
        e = eps.to_parent()
        assert e.module is ZK
        D = I.matrix.convert_to(FF(p)).hstack(*[
            (e * om).column(domain=FF(p)) for om in ZK.basis_elements()
        ])
        W = D.columnspace().convert_to(ZZ)
        H = ZK.submodule_from_matrix(W)
        factors.append(H)
    return factors


@public
def prime_decomp(p, T=None, ZK=None, dK=None, radical=None):
    r"""
    Compute the decomposition of rational prime *p* in a number field.

    Ordinarily this should be accessed through the ``primes_above()`` method
    of an :py:class:`AlgebraicField`.

    Parameters
    ==========

    p: the rational prime whose decomposition is desired.

    T: (optional) minimal, monic polynomial defining the number field $K$ in
        which to factor.

    ZK: (optional) maximal :py:class:`Order` for $K$, if already known.
        Note that at least one of *T* or *ZK* must be provided.

    dK: (optional) the discriminant of the field $K$, if already known.

    radical: (optional) :py:class:`Submodule`, giving the nilradical mod *p*
        in the integers of $K$, if already known.

    Returns
    =======

    List of :py:class:`PrimeIdeal` instances.

    Examples
    ========

    >>> from sympy import Poly, QQ
    >>> from sympy.abc import x, theta
    >>> T = Poly(x ** 3 + x ** 2 - 2 * x + 8)
    >>> K = QQ.algebraic_field((T, theta))
    >>> print(K.primes_above(2))
    [[ (2, x**2 + 1) e=1, f=1 ], [ (2, (x**2 + 3*x + 2)/2) e=1, f=1 ], [ (2, (3*x**2 + 3*x)/2) e=1, f=1 ]]

    References
    ==========

    [1] Cohen, H. *A Course in Computational Algebraic Number Theory.*
    (See Algorithm 6.2.9.)

    """
    if T is None and ZK is None:
        raise ValueError('At least one of T or ZK must be provided.')
    if T is None:
        T = ZK.T
    radicals = {}
    if dK is None or ZK is None:
        ZK, dK = round_two(T, radicals=radicals)
    dT = T.discriminant()
    f_squared = dT // dK
    if f_squared % p != 0:
        return _prime_decomp_easy_case(p, ZK)
    radical = radical or radicals.get(p) or nilradical_mod_p(ZK, p)
    stack = [radical]
    primes = []
    while stack:
        I = stack.pop()
        N, G = _prime_decomp_compute_kernel(I, p, ZK)
        if N.n == 1:
            P = _prime_decomp_maximal_ideal(I, p, ZK)
            primes.append(P)
        else:
            I1, I2 = _prime_decomp_split_ideal(I, p, N, G, ZK)
            stack.extend([I1, I2])
    return primes
