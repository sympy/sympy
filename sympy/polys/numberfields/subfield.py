r"""
Solve the "Subfield Problem" and allied problems, for algebraic number fields.

Following Cohen (see [Cohen93] Section 4.5), we can define the problem as
follows:

* **Subfield Problem:**

  Given two number fields $\mathbb{Q}(\alhpa)$, $\mathbb{Q}(\beta)$
  via the minimal polynomials for their generators $\alpha$ and $\beta$, decide
  whether one field is isomorphic to a subfield of the other.

From a solution to this problem flow solutions to the following problems as
well:

* **Primitive Element Problem:**

  Given several algebraic numbers
  $\alpha_1, \ldots, \alpha_m$, compute a single algebraic number $\theta$
  such that $\mathbb{Q}(\alpha_1, \ldots, \alpha_m) = \mathbb{Q}(\theta)$.

* **Field Isomorphism Problem:**

  Decide whether two number fields
  $\mathbb{Q}(\alpha)$, $\mathbb{Q}(\beta)$ are isomorphic.

* **Field Membership Problem:**

  Given two algebraic numbers $\alpha$,
  $\beta$, decide whether $\alpha \in \mathbb{Q}(\beta)$, and if so write
  $\alpha = f(\beta)$ for some $f(x) \in \mathbb{Q}[x]$.

This module provides solutions to all these problems.
"""

from sympy.core.add import Add
from sympy.core.numbers import AlgebraicNumber
from sympy.core.singleton import S
from sympy.core.symbol import Dummy
from sympy.core.sympify import sympify
from sympy.ntheory import sieve
from sympy.polys.densetools import dup_eval
from sympy.polys.domains import QQ
from sympy.polys.numberfields.minpoly import _choose_factor, minimal_polynomial
from sympy.polys.polyerrors import IsomorphismFailed
from sympy.polys.polytools import Poly, PurePoly, factor_list
from sympy.utilities import public

from mpmath import pslq, mp


def is_isomorphism_possible(a, b):
    """Returns `True` if there is a chance for isomorphism. """
    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if m % n != 0:
        return False

    if n == m:
        return True

    da = a.minpoly.discriminant()
    db = b.minpoly.discriminant()

    i, k, half = 1, m//n, db//2

    while True:
        p = sieve[i]
        P = p**k

        if P > half:
            break

        if ((da % p) % 2) and not (db % P):
            return False

        i += 1

    return True


def field_isomorphism_pslq(a, b):
    """Construct field isomorphism using PSLQ algorithm. """
    if not a.root.is_real or not b.root.is_real:
        raise NotImplementedError("PSLQ doesn't support complex coefficients")

    f = a.minpoly
    g = b.minpoly.replace(f.gen)

    n, m, prev = 100, b.minpoly.degree(), None

    for i in range(1, 5):
        A = a.root.evalf(n)
        B = b.root.evalf(n)

        basis = [1, B] + [ B**i for i in range(2, m) ] + [A]

        dps, mp.dps = mp.dps, n
        coeffs = pslq(basis, maxcoeff=int(1e10), maxsteps=1000)
        mp.dps = dps

        if coeffs is None:
            break

        if coeffs != prev:
            prev = coeffs
        else:
            break

        coeffs = [S(c)/coeffs[-1] for c in coeffs[:-1]]

        while not coeffs[-1]:
            coeffs.pop()

        coeffs = list(reversed(coeffs))
        h = Poly(coeffs, f.gen, domain='QQ')

        if f.compose(h).rem(g).is_zero:
            d, approx = len(coeffs) - 1, 0

            for i, coeff in enumerate(coeffs):
                approx += coeff*B**(d - i)

            if A*approx < 0:
                return [ -c for c in coeffs ]
            else:
                return coeffs
        elif f.compose(-h).rem(g).is_zero:
            return [ -c for c in coeffs ]
        else:
            n *= 2

    return None


def field_isomorphism_factor(a, b):
    """Construct field isomorphism via factorization. """
    _, factors = factor_list(a.minpoly, extension=b)

    for f, _ in factors:
        if f.degree() == 1:
            coeffs = f.rep.TC().to_sympy_list()
            d, terms = len(coeffs) - 1, []

            for i, coeff in enumerate(coeffs):
                terms.append(coeff*b.root**(d - i))

            root = Add(*terms)

            if (a.root - root).evalf(chop=True) == 0:
                return coeffs

            if (a.root + root).evalf(chop=True) == 0:
                return [-c for c in coeffs]

    return None


@public
def field_isomorphism(a, b, *, fast=True):
    """Construct an isomorphism between two number fields. """
    a, b = sympify(a), sympify(b)

    if not a.is_AlgebraicNumber:
        a = AlgebraicNumber(a)

    if not b.is_AlgebraicNumber:
        b = AlgebraicNumber(b)

    if a == b:
        return a.coeffs()

    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if n == 1:
        return [a.root]

    if m % n != 0:
        return None

    if fast:
        try:
            result = field_isomorphism_pslq(a, b)

            if result is not None:
                return result
        except NotImplementedError:
            pass

    return field_isomorphism_factor(a, b)


def _switch_domain(g, K):
    # An algebraic relation f(a, b) = 0 over Q can also be written
    # g(b) = 0 where g is in Q(a)[x] and h(a) = 0 where h is in Q(b)[x].
    # This function transforms g into h where Q(b) = K.
    frep = g.rep.inject()
    hrep = frep.eject(K, front=True)

    return g.new(hrep, g.gens[0])


def _linsolve(p):
    # Compute root of linear polynomial.
    c, d = p.rep.rep
    return -d/c


@public
def primitive_element(extension, x=None, *, ex=False, polys=False):
    """Construct a common number field for all extensions. """
    if not extension:
        raise ValueError("Cannot compute primitive element for empty extension")

    if x is not None:
        x, cls = sympify(x), Poly
    else:
        x, cls = Dummy('x'), PurePoly

    if not ex:
        gen, coeffs = extension[0], [1]
        g = minimal_polynomial(gen, x, polys=True)
        for ext in extension[1:]:
            _, factors = factor_list(g, extension=ext)
            g = _choose_factor(factors, x, gen)
            s, _, g = g.sqf_norm()
            gen += s*ext
            coeffs.append(s)

        if not polys:
            return g.as_expr(), coeffs
        else:
            return cls(g), coeffs

    gen, coeffs = extension[0], [1]
    f = minimal_polynomial(gen, x, polys=True)
    K = QQ.algebraic_field((f, gen))  # incrementally constructed field
    reps = [K.unit]  # representations of extension elements in K
    for ext in extension[1:]:
        p = minimal_polynomial(ext, x, polys=True)
        L = QQ.algebraic_field((p, ext))
        _, factors = factor_list(f, domain=L)
        f = _choose_factor(factors, x, gen)
        s, g, f = f.sqf_norm()
        gen += s*ext
        coeffs.append(s)
        K = QQ.algebraic_field((f, gen))
        h = _switch_domain(g, K)
        erep = _linsolve(h.gcd(p))  # ext as element of K
        ogen = K.unit - s*erep  # old gen as element of K
        reps = [dup_eval(_.rep, ogen, K) for _ in reps] + [erep]

    H = [_.rep for _ in reps]
    if not polys:
        return f.as_expr(), coeffs, H
    else:
        return f, coeffs, H


@public
def to_number_field(extension, theta=None, *, gen=None):
    """Express `extension` in the field generated by `theta`. """
    if hasattr(extension, '__iter__'):
        extension = list(extension)
    else:
        extension = [extension]

    if len(extension) == 1 and isinstance(extension[0], tuple):
        return AlgebraicNumber(extension[0])

    minpoly, coeffs = primitive_element(extension, gen, polys=True)
    root = sum([ coeff*ext for coeff, ext in zip(coeffs, extension) ])

    if theta is None:
        return AlgebraicNumber((minpoly, root))
    else:
        theta = sympify(theta)

        if not theta.is_AlgebraicNumber:
            theta = AlgebraicNumber(theta, gen=gen)

        coeffs = field_isomorphism(root, theta)

        if coeffs is not None:
            return AlgebraicNumber(theta, coeffs)
        else:
            raise IsomorphismFailed(
                "%s is not in a subfield of %s" % (root, theta.root))
