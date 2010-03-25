"""Tools and arithmetics for monomials of distributed polynomials. """

from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.core.basic import C

from sympy.utilities import all, any, cythonized

def monomials(variables, degree):
    """Generate a set of monomials of the given total degree or less.

       Given a set of variables `V` and a total degree `N` generate a set
       of monomials of degree at most `N`. The total number of monomials
       is defined as `(#V + N)! / (#V! N!)`, so is huge.

       For example if we would like to generate a dense polynomial of
       a total degree `N = 50` in 5 variables, assuming that exponents
       and all of coefficients are 32-bit long and stored in an array
       we would need almost 80 GiB of memory! Fortunately most
       polynomials, that we will encounter, are sparse.

       For example consider monomials in variables `x` and `y`::

           >>> from sympy import monomials
           >>> from sympy.abc import x, y

           >>> sorted(monomials([x, y], 2))
           [1, x, y, x**2, y**2, x*y]

           >>> sorted(monomials([x, y], 3))
           [1, x, y, x**2, x**3, y**2, y**3, x*y, x*y**2, y*x**2]

    """
    if not variables:
        return set([S.One])
    else:
        x, tail = variables[0], variables[1:]

        monoms = monomials(tail, degree)

        for i in range(1, degree+1):
            monoms |= set([ x**i * m for m in monomials(tail, degree-i) ])

        return monoms

def monomial_count(V, N):
    """Computes the number of monomials of degree `N` in `#V` variables.

       The number of monomials is given as `(#V + N)! / (#V! N!)`, e.g.::

           >>> from sympy import monomials, monomial_count
           >>> from sympy.abc import x, y

           >>> monomial_count(2, 2)
           6

           >>> M = monomials([x, y], 2)

           >>> sorted(M)
           [1, x, y, x**2, y**2, x*y]
           >>> len(M)
           6

    """
    return C.Factorial(V + N) / C.Factorial(V) / C.Factorial(N)

def monomial_lex_cmp(a, b):
    return cmp(a, b)

def monomial_grlex_cmp(a, b):
    return cmp(sum(a), sum(b)) or cmp(a, b)

def monomial_grevlex_cmp(a, b):
    return cmp(sum(a), sum(b)) or cmp(tuple(reversed(b)), tuple(reversed(a)))

_monomial_order = {
    'lex'     : monomial_lex_cmp,
    'grlex'   : monomial_grlex_cmp,
    'grevlex' : monomial_grevlex_cmp,
}

def monomial_cmp(order):
    """Returns a function defining admissible order on monomials.

       Currently supported orderings are:

       1. lex       - lexicographic order
       2. grlex     - graded lexicographic order
       3. grevlex   - reversed graded lexicographic order

    """
    try:
        return _monomial_order[order]
    except KeyError:
        raise ValueError("expected valid monomial order, got %s" % order)

@cythonized("a,b")
def monomial_mul(A, B):
    """Multiplication of tuples representing monomials.

       Lets multiply `x**3*y**4*z` with `x*y**2`::

           >>> from sympy.polys.monomialtools import monomial_mul

           >>> monomial_mul((3, 4, 1), (1, 2, 0))
           (4, 6, 1)

       which gives `x**4*y**5*z`.

    """
    return tuple([ a + b for a, b in zip(A, B) ])

@cythonized("a,b,c")
def monomial_div(A, B):
    """Division of tuples representing monomials.

       Lets divide `x**3*y**4*z` by `x*y**2`::

           >>> from sympy.polys.monomialtools import monomial_div

           >>> monomial_div((3, 4, 1), (1, 2, 0))
           (2, 2, 1)

       which gives `x**2*y**2*z`. However::

           >>> monomial_div((3, 4, 1), (1, 2, 2)) is None
           True

       `x*y**2*z**2` does not divide `x**3*y**4*z`.

    """
    C = [ a - b for a, b in zip(A, B) ]

    if all([ c >= 0 for c in C ]):
        return tuple(C)
    else:
        return None

@cythonized("a,b")
def monomial_gcd(A, B):
    """Greatest common divisor of tuples representing monomials.

       Lets compute GCD of `x**3*y**4*z` and `x*y**2`::

           >>> from sympy.polys.monomialtools import monomial_gcd

           >>> monomial_gcd((3, 4, 1), (1, 2, 0))
           (1, 2, 0)

       which gives `x*y**2`.

    """
    return tuple([ min(a, b) for a, b in zip(A, B) ])

@cythonized("a,b")
def monomial_lcm(A, B):
    """Least common multiple of tuples representing monomials.

       Lets compute LCM of `x**3*y**4*z` and `x*y**2`::

           >>> from sympy.polys.monomialtools import monomial_lcm

           >>> monomial_lcm((3, 4, 1), (1, 2, 0))
           (3, 4, 1)

       which gives `x**3*y**4*z`.

    """
    return tuple([ max(a, b) for a, b in zip(A, B) ])

@cythonized("i,n")
def monomial_max(*monoms):
    """Returns maximal degree for each variable in a set of monomials.

       Consider monomials `x**3*y**4*z**5`, `y**5*z` and `x**6*y**3*z**9`.
       We wish to find out what is the maximal degree for each of `x`, `y`
       and `z` variables::

           >>> from sympy.polys.monomialtools import monomial_max

           >>> monomial_max((3,4,5), (0,5,1), (6,3,9))
           (6, 5, 9)

    """
    M = list(monoms[0])

    for N in monoms[1:]:
        for i, n in enumerate(N):
            M[i] = max(M[i], n)

    return tuple(M)

@cythonized("i,n")
def monomial_min(*monoms):
    """Returns minimal degree for each variable in a set of monomials.

       Consider monomials `x**3*y**4*z**5`, `y**5*z` and `x**6*y**3*z**9`.
       We wish to find out what is the minimal degree for each of `x`, `y`
       and `z` variables::

           >>> from sympy.polys.monomialtools import monomial_min

           >>> monomial_min((3,4,5), (0,5,1), (6,3,9))
           (0, 3, 1)

    """
    M = list(monoms[0])

    for N in monoms[1:]:
        for i, n in enumerate(N):
            M[i] = min(M[i], n)

    return tuple(M)

