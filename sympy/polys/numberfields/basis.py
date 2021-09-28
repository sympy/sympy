from sympy.polys import Poly
from sympy.polys.domains import GF, ZZ
from sympy.polys.numberfields.forms import HNF, StandardRep
from sympy.polys.numberfields.modules import (
    EndomorphismRing, ModuleWithDenominator, ModuleEndomorphism, ModuleHomomorphism)
from sympy.polys.numberfields.utilities import extract_fundamental_discriminant
from sympy.utilities import public


def _apply_Dedekind_criterion(T, p):
    """
    Apply the "Dedekind criterion" to test whether the order needs to be
    enlarged relative to a given prime p.
    """
    x = T.gen
    T_bar = Poly(T, modulus=p)
    lc, fl = T_bar.factor_list()
    assert lc == 1
    g_bar = Poly(1, x, modulus=p)
    for ti_bar, _ in fl:
        g_bar *= ti_bar
    h_bar = T_bar // g_bar
    g = Poly(g_bar, domain=ZZ)
    h = Poly(h_bar, domain=ZZ)
    f = (g * h - T) // p
    f_bar = Poly(f, modulus=p)
    Z_bar = f_bar
    for b in [g_bar, h_bar]:
        Z_bar = Z_bar.gcd(b)
    U_bar = T_bar // Z_bar
    m = Z_bar.degree()
    return U_bar, m


@public
def nilradical_mod_p(H, p, q=None):
    n = H.n
    if q is None:
        q = p
        while q < n:
            q *= p
    phi = ModuleEndomorphism(H, lambda x: x**q, modulus=p)
    return phi.kernel()


def _second_enlargement(H, p, q):
    Ip = nilradical_mod_p(H, p, q=q)
    B = ModuleWithDenominator(H.W * Ip.W, H.d)
    G = B + p*H
    E = G.endomorphism_ring()
    phi = ModuleHomomorphism(H, E, lambda x: E.inner_endomorphism(x), modulus=p)
    gamma = phi.kernel()
    G = ModuleWithDenominator(H.W * gamma.W, H.d * p)
    return G + H, Ip


@public
def round_two(T, radicals=None):
    """
    Carry out Zassenhaus's "Round 2" algorithm on a monic irreducible polynomial
    T over ZZ. This computes an integral basis and the discriminant for the
    field QQ(theta), where theta is a root of T.

    Return (ZK, dK), where:
      ZK is an HNF instance representing the integral basis.
      dK is the discriminant of the field K = QQ(theta).

    Ordinarily this function need not be called directly, as one can instead
    access the :ref:`integral_basis` and :ref:`discriminant` properties of an
    :py:class:`AlgebraicField`.

    Parameters
    ----------
    T: minimal monic polynomial over ZZ defining the number field.
    radicals: optional way for any p-radicals (if computed) to be returned by
      reference. If desired, pass an empty dictionary. If the algorithm reaches
      the point where it computes the nilradical mod p of ZK, then a GF(p)-basis
      for this ideal will be stored in this dictionary under the key p. This can
      be useful for other algorithms, such as prime decomposition.

    Examples
    --------

    Calling directly:
    >>> from sympy import Poly
    >>> from sympy.abc import x
    >>> from sympy.polys.numberfields.basis import round_two
    >>> T = Poly(x ** 3 + x ** 2 - 2 * x + 8)
    >>> round_two(T)
    (HNF(Poly(x**3 + x**2 - 2*x + 8, x, domain='ZZ'), Matrix([
    [2, 0, 0],
    [0, 2, 1],
    [0, 0, 1]]), 2), -503)

    Working through an AlgebraicField:
    >>> from sympy import Poly, QQ
    >>> from sympy.abc import x, theta
    >>> T = Poly(x ** 3 + x ** 2 - 2 * x + 8)
    >>> K = QQ.algebraic_field((T, theta))
    >>> K.integral_basis
    Matrix([
    [1, 0,   0],
    [0, 1, 1/2],
    [0, 0, 1/2]])
    >>> K.discriminant
    -503
    >>> K.integral_basis_HNF
    HNF(Poly(x**3 + x**2 - 2*x + 8, x, domain='ZZ'), Matrix([
    [2, 0, 0],
    [0, 2, 1],
    [0, 0, 1]]), 2)

    References
    ==========

    [1] Cohen, H. A Course in Computational Algebraic Number Theory.

    """
    if (   not T.is_univariate
        or not T.is_irreducible
        or not T.is_monic
        or not T.domain == ZZ):
        raise ValueError('Round 2 requires a monic irreducible univariate polynomial over ZZ.')
    n = T.degree()
    D = T.discriminant()
    # D must be 0 or 1 mod 4 (see Cohen Sec 4.4), which ensures we can write
    # it in the form D = D_0 * F**2, where D_0 is 1 or a fundamental discriminant.
    _, F = extract_fundamental_discriminant(D)
    H = HNF.for_power_basis(T)
    nilrad = None
    while F:
        # Next prime:
        for p in F: break
        U_bar, m = _apply_Dedekind_criterion(T, p)
        if m == 0:
            del F[p]
            continue
        # In the first enlargement of the order spanned by the current basis, the
        # order is always ZZ[theta], and the enlargement can be done in a simpler
        # way in this case:
        U = StandardRep.from_poly(T, Poly(U_bar, domain=ZZ))
        # TODO:
        #  Theory says only first m columns of second term below are needed.
        #  Could be slightly more efficient to use only those.
        H = H + U//p*H
        if F[p] <= m:
            del F[p]
            continue
        # A second, and possibly more, enlargements for p will be needed.
        # These enlargements require a more involved procedure.
        q = p
        while q < n:
            q *= p
        H1, nilrad = _second_enlargement(H, p, q)
        while H1 != H:
            H = H1
            H1, nilrad = _second_enlargement(H, p, q)
        del F[p]
    # Note: We do not store all nilradicals mod p, only the very last. This is
    # because, unless computed against the entire integral basis, it might not
    # be accurate. (In other words, if H was not already equal to ZK when we
    # passed it to `_second_enlargement`, then we can't trust the nilradical
    # so computed.) Example: if T(x) = x ** 3 + 15 * x ** 2 - 9 * x + 13, then
    # F is divisible by 2, 3, and 7, and the nilradical mod 2 as computed above
    # will not be accurate for the full, maximal ideal ZK.
    if nilrad is not None and isinstance(radicals, dict):
        radicals[p] = nilrad
    dK = (D * H.W.det()**2) // H.d**(2*n)
    return H, dK


def _cover_more():
    """
    For testing purposes:
    Brute force search for a cubic poly that requires a second enlargement
    in Round 2. (Must set a breakpoint and run in debugger to catch it.)
    """
    from sympy.abc import x
    from sympy.polys.numberfields.utilities import coeff_search
    search = coeff_search(3, 1)
    for c in search:
        T = Poly(x ** 3 + sum(c[i] * x ** (2 - i) for i in range(3)))
        if T.is_irreducible:
            round_two(T)
