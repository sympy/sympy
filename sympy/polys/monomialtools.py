"""Tools and arithmetics for monomials of distributed polynomials. """

from sympy.core import S, C, Symbol, Mul, Tuple
from sympy.core.basic import PicklableWithSlots
from sympy.utilities import cythonized
from sympy.polys.polyerrors import ExactQuotientFailed

def monomials(variables, degree):
    r"""
    Generate a set of monomials of the given total degree or less.

    Given a set of variables `V` and a total degree `N` generate
    a set of monomials of degree at most `N`. The total number of
    monomials is huge and is given by the following formula:

    .. math::

        \frac{(\#V + N)!}{\#V! N!}

    For example if we would like to generate a dense polynomial of
    a total degree `N = 50` in 5 variables, assuming that exponents
    and all of coefficients are 32-bit long and stored in an array we
    would need almost 80 GiB of memory! Fortunately most polynomials,
    that we will encounter, are sparse.

    Examples
    ========

    Consider monomials in variables `x` and `y`::

        >>> from sympy import monomials
        >>> from sympy.abc import x, y

        >>> sorted(monomials([x, y], 2))
        [1, x, y, x**2, y**2, x*y]

        >>> sorted(monomials([x, y], 3))
        [1, x, y, x**2, x**3, y**2, y**3, x*y, x*y**2, x**2*y]

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
    r"""
    Computes the number of monomials.

    The number of monomials is given by the following formula:

    .. math::

        \frac{(\#V + N)!}{\#V! N!}

    where `N` is a total degree and `V` is a set of variables.

    Examples
    ========

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
    return C.factorial(V + N) / C.factorial(V) / C.factorial(N)

class MonomialOrder(object):
    """Base class for monomial orderings. """

    alias = None

    def key(self, monomial):
        raise NotImplementedError

    def __str__(self):
        return self.alias

    def __call__(self, monomial):
        return self.key(monomial)

class LexOrder(MonomialOrder):
    """Lexicographic order of monomials. """

    alias = 'lex'

    def key(self, monomial):
        return monomial

class GradedLexOrder(MonomialOrder):
    """Graded lexicographic order of monomials. """

    alias = 'grlex'

    def key(self, monomial):
        return (sum(monomial), monomial)

class ReversedGradedLexOrder(MonomialOrder):
    """Reversed graded lexicographic order of monomials. """

    alias = 'grevlex'

    def key(self, monomial):
        return (sum(monomial), tuple(reversed([-m for m in monomial])))

lex = LexOrder()
grlex = GradedLexOrder()
grevlex = ReversedGradedLexOrder()

_monomial_key = {
    'lex'     : lex,
    'grlex'   : grlex,
    'grevlex' : grevlex,
}

def monomial_key(order=None):
    """
    Return a function defining admissible order on monomials.

    The result of a call to :func:`monomial_key` is a function which should
    be used as a key to :func:`sorted` built-in function, to provide order
    in a set of monomials of the same length.

    Currently supported monomial orderings are:

    1. lex       - lexicographic order (default)
    2. grlex     - graded lexicographic order
    3. grevlex   - reversed graded lexicographic order

    If the input argument is not a string but has ``__call__`` attribute,
    then it will pass through with an assumption that the callable object
    defines an admissible order on monomials.

    """
    if order is None:
        return lex

    if isinstance(order, Symbol):
        order = str(order)

    if isinstance(order, str):
        try:
            return _monomial_key[order]
        except KeyError:
            raise ValueError("supported monomial orderings are 'lex', 'grlex' and 'grevlex', got %r" % order)
    elif hasattr(order, '__call__'):
        return order
    else:
        raise ValueError("monomial ordering specification must be a string or a callable, got %s" % order)

@cythonized("a,b")
def monomial_mul(A, B):
    """
    Multiplication of tuples representing monomials.

    Lets multiply `x**3*y**4*z` with `x*y**2`::

        >>> from sympy.polys.monomialtools import monomial_mul

        >>> monomial_mul((3, 4, 1), (1, 2, 0))
        (4, 6, 1)

    which gives `x**4*y**5*z`.

    """
    return tuple([ a + b for a, b in zip(A, B) ])

@cythonized("a,b,c")
def monomial_div(A, B):
    """
    Division of tuples representing monomials.

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

    if all(c >= 0 for c in C):
        return tuple(C)
    else:
        return None

@cythonized("a,b")
def monomial_gcd(A, B):
    """
    Greatest common divisor of tuples representing monomials.

    Lets compute GCD of `x*y**4*z` and `x**3*y**2`::

        >>> from sympy.polys.monomialtools import monomial_gcd

        >>> monomial_gcd((1, 4, 1), (3, 2, 0))
        (1, 2, 0)

    which gives `x*y**2`.

    """
    return tuple([ min(a, b) for a, b in zip(A, B) ])

@cythonized("a,b")
def monomial_lcm(A, B):
    """
    Least common multiple of tuples representing monomials.

    Lets compute LCM of `x*y**4*z` and `x**3*y**2`::

        >>> from sympy.polys.monomialtools import monomial_lcm

        >>> monomial_lcm((1, 4, 1), (3, 2, 0))
        (3, 4, 1)

    which gives `x**3*y**4*z`.

    """
    return tuple([ max(a, b) for a, b in zip(A, B) ])

@cythonized("i,n")
def monomial_max(*monoms):
    """
    Returns maximal degree for each variable in a set of monomials.

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
    """
    Returns minimal degree for each variable in a set of monomials.

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

class Monomial(PicklableWithSlots):
    """Class representing a monomial, i.e. a product of powers. """

    __slots__ = ['exponents', 'gens']

    def __init__(self, exponents, gens=None):
        self.exponents = tuple(exponents)
        self.gens = gens

    def rebuild(self, exponents, gens=None):
        return self.__class__(exponents, gens or self.gens)

    def __len__(self):
        return len(self.exponents)

    def __iter__(self):
        return iter(self.exponents)

    def __getitem__(self, item):
        return self.exponents[item]

    def __hash__(self):
        return hash((self.__class__.__name__, self.exponents, self.gens))

    def __str__(self):
        if self.gens:
            return "*".join([ "%s**%s" % (gen, exp) for gen, exp in zip(self.gens, self.exponents) ])
        else:
            return "%s(%s)" % (self.__class__.__name__, self.exponents)

    def as_expr(self, *gens):
        """Convert a monomial instance to a SymPy expression. """
        gens = gens or self.gens

        if not gens:
            raise ValueError("can't convert %s to an expression without generators" % self)

        return Mul(*[ gen**exp for gen, exp in zip(gens, self.exponents) ])

    def __eq__(self, other):
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return False

        return self.exponents == exponents

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, other):
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return NotImplementedError

        return self.rebuild(monomial_mul(self.exponents, exponents))

    def __div__(self, other):
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return NotImplementedError

        result = monomial_div(self.exponents, exponents)

        if result is not None:
            return self.rebuild(result)
        else:
            raise ExactQuotientFailed(self, Monomial(other))

    __floordiv__ = __truediv__ = __div__

    def __pow__(self, other):
        n = int(other)

        if not n:
            return self.rebuild([0]*len(self))
        elif n > 0:
            exponents = self.exponents

            for i in xrange(1, n):
                exponents = monomial_mul(exponents, self.exponents)

            return self.rebuild(exponents)
        else:
            raise ValueError("a non-negative integer expected, got %s" % other)

    def gcd(self, other):
        """Greatest common divisor of monomials. """
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            raise TypeError("an instance of Monomial class expected, got %s" % other)

        return self.rebuild(monomial_gcd(self.exponents, exponents))

    def lcm(self, other):
        """Least common multiple of monomials. """
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            raise TypeError("an instance of Monomial class expected, got %s" % other)

        return self.rebuild(monomial_lcm(self.exponents, exponents))
