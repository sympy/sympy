from sympy.polys.monomials import monomial_ngcd
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.external.gmpy import MPZ
    from sympy.polys.sparsetools import smp
    from sympy.polys.domains.domain import Domain, Er

def gcd_terms(polynomials: list[smp[Er]], domain: Domain[Er]) -> smp[Er]:
    """
    Returns the greatest common divisor (GCD) of all terms in a list of
    polynomials p with respect to a given ring and domain.

    Examples
    ========

    >>> from sympy.polys.rings import gcd_terms
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)

    >>> p1 = 2*x**3*y**2 - 2*x**2*y**3
    >>> p2 = 4*x**2*y**3 - 4*x**3*y**2
    >>> polynomials = [p1, p2]
    >>> ring = p1.ring
    >>> domain = ring.domain
    >>> R(gcd_terms(polynomials, domain))
    2*x**2*y**2
    >>> p1.gcd(p2) # Shows the difference between the gcd_terms and gcd
    2*x**3*y**2 - 2*x**2*y**3

    """
    monomials = set()
    coeffs = set()

    for pi in polynomials:
        for monomial, coeff in pi.items():
            monomials.add(monomial)
            coeffs.add(coeff)

    monom_gcd = monomial_ngcd(list(monomials))
    coeff_gcd = domain.gcdn(coeffs)

    return {monom_gcd: coeff_gcd}


def gcd_prs(f, g):
    """
    Returns the greatest common divisor (GCD) of two polynomials using the
    Polynomial Resultant Sequences (PRS) method.

    Examples
    ========

    >>> from sympy.polys.rings import gcd_prs
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)

    >>> f = 4*x**2 + 8*x + 10*y
    >>> g = 2*x**3 + 4*x**2*y + 2*x*y**2
    >>> gcd_prs(f, g)
    2

    """
    K = f.ring.domain

    x = f.main_variable()

    c1, pp1 = f.primitive_wrt(x)
    c2, pp2 = g.primitive_wrt(x)

    c = gcd_extract([c1, c2],"prs")

    h = pp1.subresultants(pp2, x)[-1]
    _, h = h.primitive_wrt(x)

    if not K.is_Field:
        c *= K.canonical_unit(h.LC)

    h = c * h

    if K.is_Field:
        h = h.monic()

    return h


def gcd_extract(polynomials, tag):
    """
    Returns the greatest common divisor (GCD) of a list of polynomials.

    Examples
    ========

    >>> from sympy.polys.rings import gcd_extract
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)

    >>> polynomials = x - y, x - y
    >>> gcd_extract(polynomials, "heugcd")
    x - y

    """
    ring = polynomials[0].ring
    domain = ring.domain

    if any(len(pi) == 1 for pi in polynomials):
        return gcd_terms(polynomials, domain)

    polynomials, monom_gcd = monomial_extract(polynomials)

    polynomials, common_symbols = _gcd_preprocess_polys(polynomials)

    if tag == "prs":
        if len(polynomials) == 1:
            gcd = polynomials[0]
            if monom_gcd is not None:
                gcd = gcd * monom_gcd

        gcd = polynomials[0]
        for pi in polynomials[1:]:
            gcd = gcd_prs(gcd, pi)
            if gcd == domain.one:
                break

        if monom_gcd is not None:
            gcd = gcd * monom_gcd

        if domain.is_Field:
            gcd = gcd.monic()
        else:
            gcd *= domain.canonical_unit(gcd.LC)

        return gcd

    elif tag == "heugcd":
        if len(polynomials) == 1:
            gcd = polynomials[0]
            if monom_gcd is not None:
                gcd = gcd * monom_gcd

            return gcd

        gcd = polynomials[0]
        for pi in polynomials[1:]:
            gcd = heugcd(gcd, pi)[0]
            if gcd == domain.one:
                break

        if monom_gcd is not None:
            gcd = gcd * monom_gcd

        return gcd


def _gcd_preprocess_polys(polynomials, dom):
    """
    Simplify a list of polynomials whose gcd is wanted. Returns a possibly
    longer list of simpler polynomials having the same gcd as the input.

    This is done by eliminating symbols that can not be part of the gcd
    because they do not appear in each item of the input. In the output
    list, all items have exactly the same symbols. The set of those symbols
    is also returned.

    Examples
    ========

    >>> from sympy.core import ordered
    >>> from sympy.polys.rings import _gcd_preprocess_polys
    >>> from sympy import ring, ZZ
    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2 - y**2
    >>> g = x**2 - 2*x*y + y**2
    >>> h = x - y
    >>> polynomials = [f, g, h]
    >>> result = _gcd_preprocess_polys(polynomials)
    >>> list(ordered(result[0])), result[1] # Set ordering is random
    ([x - y, x**2 - y**2, x**2 - 2*x*y + y**2], {0, 1})

    """
    all_polys = polynomials

    while True:
        # Quick exits are most efficient if we start from the simplest polys
        polynomials = sorted(set(all_polys), key=len)

        # aggiungi euristica che controlla se i pol hanno lunghezza 1 come fatto in fondo

        # Find the intersection of symbols for each poly:
        common = polynomials[0].free_variables()
        allsame = True
        for pi in polynomials[1:]:

            # Quick exit
            if not common:
                ring = polynomials[0].ring
                domain = ring.domain
                gcd = gcd_terms(polynomials, dom) #I think that here the gcd should be scalar, and we shouldn't call gcd_terms
                return [gcd], None

            syms = pi.free_variables()
            if allsame and syms != common:
                allsame = False
            common &= syms

        # The loop is complete if they all have the same symbols.
        if allsame:
            return polynomials, common

        # Extract coefficients as polys containing only the common symbols.
        all_polys = []
        for i, pi in enumerate(polynomials):
            coeffs_i = pi.coeff_split(pi.free_variables() - common)
            all_polys.extend(coeffs_i)

            # Quick exit:
            if any(len(c) == 1 for c in coeffs_i):
                ring = polynomials[0].ring
                domain = ring.domain
                gcd = gcd_terms((all_polys + polynomials[i+1:]), dom)
                return [gcd], None


def _gcd_ring(self, g):
    """Helper function for ``_gcd`` method."""

    f = self
    K = f.ring.domain

    if not K.is_Exact:
        K_exact = K.get_exact()

        ring_approx = f.ring
        ring_exact = K_exact[ring_approx.symbols].ring

        f = f.set_ring(ring_exact)
        g = g.set_ring(ring_exact)

        h, cff, cfg = f._gcd_ring(g)

        h = h.set_ring(ring_approx)
        cff = cff.set_ring(ring_approx)
        cfg = cfg.set_ring(ring_approx)

        return h, cff, cfg
    elif K.is_Field:
        if K.is_QQ:
            return f._gcd_QQ(g)

        gcd = gcd_extract((self,g), "prs")
        cff, cfg = f.quo(gcd), g.quo(gcd)
        return gcd, cff, cfg

    else:
        if K.is_ZZ:
            return f._gcd_ZZ(g)

        gcd = gcd_extract((self,g), "prs")
        cff, cfg = f.quo(gcd), g.quo(gcd)
        return gcd, cff, cfg


def primitive_wrt(self, x):
    """
    Returns the content and primitive part of a polynomial with respect to a
    specified variable.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x, y = ring("x, y", ZZ)
    >>> p = 6*x**2*y - 9*x**3*y**2 + 3*x*y
    >>> p.primitive_wrt(x)
    (3*y, -3*x**3*y + 2*x**2 + x)

    >>> p.primitive() # Distinguishing primitive and primitive_wrt outcomes
    (3, -3*x**3*y**2 + 2*x**2*y + x*y)
    """
    p = self
    coeffs = p.coeff_split({x})
    cont = gcd_extract(coeffs, "prs")
    prim = p.exquo(cont)
    return cont, prim


def zippel_preprocess(A, B):
    """
    checks if in any variable the gcd is monic
    if it doesn't find any variable with monic gcd
    it extracts content wrt the first variable
    and proceeds to compute the gcd
    maybe it is better to call the generic sympy gcd rather than
    modgcd_multivariate once the whole dispatching heuristic is implemented
    """
    ring = A.ring
    ngens = ring.ngens
    syms = ring.symbols
    t = None

    for i in range(ngens):
        m_A = max(el[i]for el in A.keys())
        m_B = max(el[i] for el in B.keys())
        if t == None and (m_A != 0 or m_B != 0):
            t = i
        lc_A = A.coeff_wrt(i, m_A)
        lc_B = B.coeff_wrt(i, m_B)
        if len(lc_A) == 1 or len(lc_B) == 1:
            if i == 0:
                gcd, _, _ = modgcd_multivariate(A, B, "zippel")
                return gcd
            new_syms = (syms[i],) + syms[:i] + syms[i+1:]
            new_ring = ring.clone(symbols = new_syms)
            A_ = A.set_ring(new_ring)
            B_ = B.set_ring(new_ring)
            gcd, _, _ = modgcd_multivariate(A_, B_, "zippel")
            return gcd.set_ring(ring)
    if t == None:
        gcd, _, _ = modgcd_multivariate(A, B, "zippel")
        return gcd

    A_cont, A_prim = A.primitive_wrt(t)
    B_cont, B_prim = B.primitive_wrt(t)
    cont_gcd = zippel_preprocess(A_cont, B_cont)
    gcd,_ , _ = modgcd_multivariate(A_prim, B_prim, "zippel")
    return gcd * cont_gcd


def gcd_extract(v):

    # faster than computing the gcd iteratively
    # when gcd is computed with Zippel

    if len(v) == 1:
        return v[1]
    if len(v) == 2:
        return zippel_preprocess(v_sor[0], v_sor[1])
    v_sor = sorted(v, key=lambda p: len(p))
    A = v_sor[1]
    while True:
        for el in v_sor[2:]:
            k = random.randint(1, 10009)
            A += el * k
        gcd = zippel_preprocess(v_sor[0], A)
        restart = False
        for pol in v:
            if pol.rem(gcd) != 0:
                restart = True
                break
        if not restart:
            break

    return gcd
