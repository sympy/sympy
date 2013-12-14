from __future__ import print_function, division

from sympy.polys.polytools import factor_list


def dispersionset(p, q=None, gen=None):
    r"""Compute the 'dispersion set' of two polynomials.

    For two polynomials `f(x)` and `g(x)` with `deg f > 0`
    and `deg g > 0` the dispersion set `J(f,g)` is defined as:

    .. math::
        J(f, g) & := \{a \in Z^{+} | gcd(f(x), g(x+a)) \neq 1\} \\
                &  = \{a \in Z^{+} | deg gcd(f(x), g(x+a)) \geq 1\}

    For a single polynomial one defines `J(f) := J(f, f)`.

    See Also
    ========

    dispersion

    References
    ==========

    ..[1]: "On the Summation of Rational Functions"
    ..[2]: "On Computing Closed Forms for Indefinite Summations"
    ..[3]: "Hypergeometric Summation: An Algorithmic Approach to Summationand Special Function Identities.
    ..[4]: "Fast Polynomial Dispersion Computation and its Application to Indefinite Summation"
    """
    # Check for valid input
    if q is None:
        q = p
    if not p.is_univariate or not q.is_univariate:
        raise ValueError("Polynomials need to be univariate")
    # We define the dispersion of constant polynomials to be zero
    if p.degree() < 1 or q.degree() < 1:
        return set([0])
    # The generator
    if gen is None:
        gen = p.gens[0]
    # Factor p and q over the rationals
    fp = factor_list(p)
    fq = factor_list(q)
    J = set([])
    # Iterate over all pairs of factors
    for s, unused in fp[1]:
        for t, unused in fq[1]:
            m = s.degree()
            n = t.degree()
            if not m == n:
                D = []
            else:
                a = s.coeff_monomial(gen**n)
                b = s.coeff_monomial(gen**(n-1))
                c = t.coeff_monomial(gen**n)
                d = t.coeff_monomial(gen**(n-1))
                j = (b*c-a*d) / (a*c*n)
                if j < 0:
                    D = []
                else:
                    if (c*s - a*t.shift(j)).is_zero:
                        D = [j]
                    else:
                        D = []
            J = J.union(D)
    return J


def dispersion(p, q=None, gen=None):
    r"""Compute the 'dispersion' of a polynomial.

    For two polynomials `f(x)` and `g(x)` with `deg f > 0`
    and `deg g > 0` the dispersion `dis(f,g)` is defined as:

    .. math::
        dis(f, g) & := max\{ J(f,g) \union \{0\}\} \\
                  &  = max\{ \{a \in Z^{+} | gcd(f(x), g(x+a)) \neq 1\} \union \{0\}\}

    and for a single polynomial `dis(f) := dis(f, f)`.

    See Also
    ========

    dispersionset

    References
    ==========

    ..[1]: "On the Summation of Rational Functions"
    ..[2]: "On Computing Closed Forms for Indefinite Summations"
    ..[3]: "Hypergeometric Summation: An Algorithmic Approach to Summationand Special Function Identities.
    ..[4]: "Fast Polynomial Dispersion Computation and its Application to Indefinite Summation"
    """
    J = dispersionset(p, q, gen)
    return max(J)
