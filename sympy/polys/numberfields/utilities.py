"""Utilities for algebraic number theory. """

from sympy.core import Basic
from sympy.utilities import public

def is_int(a):
    return isinstance(a, int) or (isinstance(a, Basic) and a.is_integer)

def is_rat(a):
    return is_int(a) or (isinstance(a, Basic) and a.is_rational)

def is_zero(a):
    return (is_int(a) or is_rat(a)) and a == 0

@public
def extract_fundamental_discriminant(a):
    """
    Given any rational integer `a` that is 0 or 1 mod 4, write a = d*f**2,
    where d is either 1 or a fundamental discriminant, and return a pair
    of dictionaries (D, F) giving the prime factorizations of d and f resp.,
    in the same format returned by `sympy.ntheory.factor_.factorint`.

    A fundamental discriminant d is != 1 and is either 1 mod 4 and squarefree,
    or is 0 mod 4 and such that d/4 is squarefree and 2 or 3 mod 4. (This is
    the same as being the discriminant of a quadratic field.)

    Examples
    ========

    >>> from sympy.polys.numberfields.utilities import extract_fundamental_discriminant
    >>> extract_fundamental_discriminant(-432)
    ({3: 1, -1: 1}, {2: 2, 3: 1})

    For comparison:
    >>> from sympy import factorint
    >>> factorint(-432)
    {2: 4, 3: 3, -1: 1}
    """
    from sympy.ntheory.factor_ import factorint
    if a % 4 not in [0, 1]:
        raise ValueError('To extract fundamental discriminant, number must be 0 or 1 mod 4.')
    if a == 0:
        return {}, {0: 1}
    if a == 1:
        return {}, {}
    a_factors = factorint(a)
    D = {}
    F = {}
    # First pass: just make d squarefree, and a/d a perfect square.
    # We'll count primes (and units! i.e. -1) that are 3 mod 4 dividing d.
    num_3_mod_4 = 0
    for p, e in a_factors.items():
        if e % 2 == 1:
            D[p] = 1
            if p % 4 == 3:
                num_3_mod_4 += 1
            if e >= 3:
                F[p] = (e - 1) // 2
        else:
            F[p] = e // 2
    # Second pass: if d is cong. to 2 or 3 mod 4, then we must steal away
    # another factor of 4 from f**2 and give it to d.
    even = 2 in D
    if even or num_3_mod_4 % 2 == 1:
        e2 = F[2]
        assert e2 > 0
        if e2 == 1:
            del F[2]
        else:
            F[2] = e2 - 1
        D[2] = 3 if even else 2
    return D, F


@public
class AlgIntPowers:
    """
    Given an algebraic integer theta by its monic minimal polynomial T over ZZ,
    this class computes representations of arbitrarily high powers of theta, as
    ZZ-linear combinations over {1, theta, ..., theta**(n-1)}, where n = deg(T).

    Optionally, the representations may be reduced w.r.t. a modulus.

    Examples
    ========

    >>> from sympy import Poly, cyclotomic_poly
    >>> from sympy.polys.numberfields.utilities import AlgIntPowers
    >>> T = Poly(cyclotomic_poly(5))
    >>> zeta_pow = AlgIntPowers(T)
    >>> zeta_pow[0]
    [1, 0, 0, 0]
    >>> zeta_pow[1]
    [0, 1, 0, 0]
    >>> zeta_pow[4]
    [-1, -1, -1, -1]
    >>> zeta_pow[24]
    [-1, -1, -1, -1]

    """

    def __init__(self, T, modulus=None):
        """
        Parameters
        ----------
        T: the monic minimal polynomial over ZZ defining the algebraic integer.
        modulus: if not None, all representations will be reduced w.r.t. this.
        """
        self.T = T
        self.modulus = modulus
        self.n = T.degree()
        self.powers_n_and_up = [[-c % self for c in reversed(T.rep.rep)][:-1]]
        self.max_so_far = self.n

    def red(self, exp):
        return exp if self.modulus is None else exp % self.modulus

    def __rmod__(self, other):
        return self.red(other)

    def compute_up_through(self, e):
        m = self.max_so_far
        if e <= m: return
        n = self.n
        r = self.powers_n_and_up
        c = r[0]
        for k in range(m+1, e+1):
            b = r[k-1-n][n-1]
            r.append(
                [c[0]*b % self] + [
                    (r[k-1-n][i-1] + c[i]*b) % self for i in range(1, n)
                ]
            )
        self.max_so_far = e

    def get(self, e):
        n = self.n
        if e < 0:
            raise ValueError('Exponent must be non-negative.')
        elif e < n:
            return [1 if i == e else 0 for i in range(n)]
        else:
            self.compute_up_through(e)
            return self.powers_n_and_up[e - n]

    def __getitem__(self, item):
        return self.get(item)


@public
def coeff_search(m, R):
    """
    Generate coefficients for searching through polynomials.

    Parameters
    ----------
    m: length of coeff list
    R: initial max abs val for coeffs (will increase later)

    Returns
    -------
    Infinite generator of lists of coefficients.

    """
    R0 = R
    c = [R] * m
    while True:
        if R == R0 or R in c or -R in c:
            yield c[:]
        j = m - 1
        while c[j] == -R:
            j -= 1
        c[j] -= 1
        for i in range(j + 1, m):
            c[i] = R
        for j in range(m):
            if c[j] != 0:
                break
        else:
            R += 1
            c = [R] * m
