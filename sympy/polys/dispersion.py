from __future__ import print_function, division

from sympy.core import S
from sympy.polys.polytools import factor_list


def dispersionset(p, q=None):
    r"""Compute the *dispersion set* of two polynomials.

    For two polynomials `f(x)` and `g(x)` with `\deg f > 0`
    and `\deg g > 0` the dispersion set `\operatorname{J}(f, g)` is defined as:

    .. math::
        \operatorname{J}(f, g)
        & := \{a \in \mathbb{N}_0 | \gcd(f(x), g(x+a)) \neq 1\} \\
        &  = \{a \in \mathbb{N}_0 | \deg \gcd(f(x), g(x+a)) \geq 1\}

    For a single polynomial one defines `\operatorname{J}(f) := \operatorname{J}(f, f)`.

    See Also
    ========

    dispersion

    References
    ==========

    1. [ManWright94]_
    2. [Koepf98]_
    3. [Abramov71]_
    4. [Man93]_
    """
    # Check for valid input
    same = False if q is not None else True
    if same:
        q = p

    if not p.is_univariate or not q.is_univariate:
        raise ValueError("Polynomials need to be univariate")

    # The generator
    if not p.gen == q.gen:
        raise ValueError("Polynomials must have the same generator")
    gen = p.gen

    # We define the dispersion of constant polynomials to be zero
    if p.degree() < 1 or q.degree() < 1:
        return set([0])

    # Factor p and q over the rationals
    fp = p.factor_list()
    fq = q.factor_list() if not same else fp

    # Iterate over all pairs of factors
    J = set([])
    for s, unused in fp[1]:
        for t, unused in fq[1]:
            m = s.degree()
            n = t.degree()
            if n != m:
                continue
            an = s.LC()
            bn = t.LC()
            if not (an - bn).is_zero:
                continue
            # Note that the roles of `s` and `t` below are switched
            # w.r.t. the original paper. This is for consistency
            # with the description in the book of W. Koepf.
            anm1 = s.coeff_monomial(gen**(m-1))
            bnm1 = t.coeff_monomial(gen**(n-1))
            alpha = (anm1 - bnm1) / S(n*bn)
            if not alpha.is_integer:
                continue
            if alpha < 0 or alpha in J:
                continue
            if n > 1 and not (s - t.shift(alpha)).is_zero:
                continue
            J.add(alpha)

    return J


def dispersion(p, q=None):
    r"""Compute the *dispersion* of polynomials.

    For two polynomials `f(x)` and `g(x)` with `\deg f > 0`
    and `\deg g > 0` the dispersion `\operatorname{dis}(f, g)` is defined as:

    .. math::
        \operatorname{dis}(f, g)
        & := \max\{ J(f,g) \cup \{0\} \} \\
        &  = \max\{ \{a \in \mathbb{N} | \gcd(f(x), g(x+a)) \neq 1\} \cup \{0\} \}

    and for a single polynomial `\operatorname{dis}(f) := \operatorname{dis}(f, f)`.

    See Also
    ========

    dispersionset

    References
    ==========

    1. [ManWright94]_
    2. [Koepf98]_
    3. [Abramov71]_
    4. [Man93]_
    """
    J = dispersionset(p, q)
    if len(J) == 0:
        # Definition for maximum of empty set
        j = S.NegativeInfinity
    else:
        j = max(J)
    return j
