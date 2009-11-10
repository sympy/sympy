
from sympy.core.mul import Mul
from sympy.core.basic import S
from sympy.core.cache import cacheit

from sympy.functions import factorial
from sympy.utilities import all, any

def monomials(variables, degree):
    """Generate a set of monomials of the given total degree or less.

       Given a set of variables V and a total degree  N generate a set
       of monomials of degree at most N. The total number of monomials
       is defined as (#V + N)! / (#V! N!) so is huge.

       For example if we would like to generate a dense polynomial of
       a total degree N = 50 in 5 variables,  assuming that exponents
       and all of coefficients are 32-bit long and stored in an array
       we would need almost 80 GiB of memory! Fortunately most
       polynomials, that we will encounter, are sparse.

       >>> from sympy.polys.monomial import monomials
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

@cacheit
def monomial_count(V, N):
    """Computes the number of monomials of degree N in #V variables.

       The number of monomials is given as (#V + N)! / (#V! N!), eg:

       >>> from sympy.polys.monomial import monomial_count, monomials
       >>> from sympy.abc import x, y

       >>> monomial_count(2, 2)
       6

       >>> M = monomials([x, y], 2)

       >>> print M
       set([1, x, y, x*y, x**2, y**2])
       >>> len(M)
       6

    """
    return factorial(V + N) / factorial(V) / factorial(N)

@cacheit
def monomial_lex_cmp(a, b):
    return cmp(a, b)

@cacheit
def monomial_grlex_cmp(a, b):
    return cmp(sum(a), sum(b)) or cmp(a, b)

@cacheit
def monomial_grevlex_cmp(a, b):
    return cmp(sum(a), sum(b)) or cmp(tuple(reversed(b)), tuple(reversed(a)))

@cacheit
def monomial_1_el_cmp(a, b):
    return cmp(a[0], b[0]) or cmp(sum(a[2:]),
        sum(b[2:])) or cmp(tuple(reversed(b[2:])), tuple(reversed(a[2:])))

_monomial_order = {
    'lex'     : monomial_lex_cmp,
    'grlex'   : monomial_grlex_cmp,
    'grevlex' : monomial_grevlex_cmp,
    '1-el'    : monomial_1_el_cmp,
}

def monomial_cmp(order):
    """Returns a function defining admissible order on monomials.

       Currently supported orderings are:

           [1] lex       -> lexicographic order
           [2] grlex     -> graded lex order
           [3] grevlex   -> reversed grlex order
           [4] 1-el      -> first elimination order

    """
    try:
        return _monomial_order[order]
    except KeyError:
        raise ValueError("Unknown monomial order: %s" % order)

@cacheit
def monomial_mul(a, b):
    """Multiplication of tuples representing monomials.

       Lets multiply x**3*y**4*z with x*y**2:

       >>> from sympy.polys.monomial import monomial_mul
       >>> monomial_mul((3, 4, 1), (1, 2, 0))
       (4, 6, 1)

       which gives x**4*y**6*z.

    """
    return tuple([ x + y for x, y in zip(a, b) ])

@cacheit
def monomial_div(a, b):
    """Division of tuples representing monomials.

       Lets divide x**3*y**4*z by x*y**2:

       >>> from sympy.polys.monomial import monomial_div
       >>> monomial_div((3, 4, 1), (1, 2, 0))
       (2, 2, 1)

       which gives x**2*y**2*z. However, nothing is obtained for the
       following,

       >>> monomial_div((3, 4, 1), (1, 2, 2))
       <BLANKLINE>

       since x*y**2*z**2 does not divide x**3*y**4*z.
    """
    result = [ x - y for x, y in zip(a, b) ]

    if all(e >= 0 for e in result):
        return tuple(result)
    else:
        return None

@cacheit
def monomial_gcd(a, b):
    """Greatest common divisor of tuples representing monomials.

       Lets compute GCD of x**3*y**4*z and x*y**2:

       >>> from sympy.polys.monomial import monomial_gcd
       >>> monomial_gcd((3, 4, 1), (1, 2, 0))
       (1, 2, 0)

       which gives x*y**2.

    """
    return tuple([ min(x, y) for x, y in zip(a, b) ])

@cacheit
def monomial_lcm(a, b):
    """Least common multiple of tuples representing monomials.

       Lets compute LCM of x**3*y**4*z and x*y**2:

       >>> from sympy.polys.monomial import monomial_lcm
       >>> monomial_lcm((3, 4, 1), (1, 2, 0))
       (3, 4, 1)

       which gives x**3*y**4*z.

    """
    return tuple([ max(x, y) for x, y in zip(a, b) ])

@cacheit
def monomial_max(*monoms):
    """Returns maximal degree for each variable in a set of monomials.

       Consider monomials x**3*y**4*z**5, y**5*z and x**6*y**3*z**9. We
       wish to find out what is the maximal degree for each of x, y, z
       variables:

       >>> from sympy.polys.monomial import monomial_max
       >>> monomial_max((3,4,5), (0,5,1), (6,3,9))
       (6, 5, 9)

    """
    return tuple(map(lambda *row: max(row), *monoms))

@cacheit
def monomial_min(*monoms):
    """Returns minimal degree for each variable in a set of monomials.

       Consider monomials x**3*y**4*z**5, y**5*z and x**6*y**3*z**9. We
       wish to find out what is the maximal degree for each of x, y, z
       variables:

       >>> from sympy.polys.monomial import monomial_min
       >>> monomial_min((3,4,5), (0,5,1), (6,3,9))
       (0, 3, 1)

    """
    return tuple(map(lambda *row: min(row), *monoms))

@cacheit
def monomial_as_basic(monom, *syms):
    """Converts tuple representing monomial to a valid sympy expression.

       Given a monomial and a list of symbols, both tuples must be of
       the same length, returns a sympy expression representing this
       monomial, eg. consider monomial (3, 2, 1) over (x, y, z):

       >>> from sympy.polys.monomial import monomial_as_basic
       >>> from sympy.abc import x, y, z

       >>> monomial_as_basic((3, 2, 1), x, y, z)
       z*x**3*y**2

    """
    return Mul(*[ b**e for b, e in zip(syms, monom) ])
