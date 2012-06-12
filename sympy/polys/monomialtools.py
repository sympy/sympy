"""Tools and arithmetics for monomials of distributed polynomials. """

from sympy.core import S, C, Symbol, Mul, Tuple
from sympy.polys.polyutils import PicklableWithSlots
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
    is_global = None

    def key(self, monomial):
        raise NotImplementedError

    def __str__(self):
        return self.alias

    def __call__(self, monomial):
        return self.key(monomial)

    def __eq__(self, other):
        return self.__class__ == other.__class__

    def __hash__(self):
        return hash(self.__class__)

    def __ne__(self, other):
        return not (self == other)

class LexOrder(MonomialOrder):
    """Lexicographic order of monomials. """

    alias = 'lex'
    is_global = True

    def key(self, monomial):
        return monomial

class GradedLexOrder(MonomialOrder):
    """Graded lexicographic order of monomials. """

    alias = 'grlex'
    is_global = True

    def key(self, monomial):
        return (sum(monomial), monomial)

class ReversedGradedLexOrder(MonomialOrder):
    """Reversed graded lexicographic order of monomials. """

    alias = 'grevlex'
    is_global = True

    def key(self, monomial):
        return (sum(monomial), tuple(reversed([-m for m in monomial])))

class ProductOrder(MonomialOrder):
    """
    A product order built from other monomial orders.

    Given (not necessarily total) orders O1, O2, ..., On, their product order
    P is defined as M1 > M2 iff there exists i such that O1(M1) = O2(M2),
    ..., Oi(M1) = Oi(M2), O{i+1}(M1) > O{i+1}(M2).

    Product orders are typically built from monomial orders on different sets
    of variables.

    ProductOrder is constructed by passing a list of pairs
    [(O1, L1), (O2, L2), ...] where Oi are MonomialOrders and Li are callables.
    Upon comparison, the Li are passed the total monomial, and should filter
    out the part of the monomial to pass to Oi.

    Examples
    ========

    We can use a lexicographic order on x_1, x_2 and also on
    y_1, y_2, y_3, and their product on {x_i, y_i} as follows:

    >>> from sympy.polys.monomialtools import lex, grlex, ProductOrder
    >>> P = ProductOrder(
    ...     (lex, lambda m: m[:2]), # lex order on x_1 and x_2 of monomial
    ...     (grlex, lambda m: m[2:]) # grlex on y_1, y_2, y_3
    ... )
    >>> P((2, 1, 1, 0, 0)) > P((1, 10, 0, 2, 0))
    True

    Here the exponent `2` of `x_1` in the first monomial
    (`x_1^2 x_2 y_1`) is bigger than the exponent `1` of `x_1` in the
    second monomial (`x_1 x_2^10 y_2^2`), so the first monomial is greater
    in the product ordering.

    >>> P((2, 1, 1, 0, 0)) < P((2, 1, 0, 2, 0))
    True

    Here the exponents of `x_1` and `x_2` agree, so the grlex order on
    `y_1, y_2, y_3` is used to decide the ordering. In this case the monomial
    `y_2^2` is ordered larger than `y_1`, since for the grlex order the degree
    of the monomial is most important.
    """

    def __init__(self, *args):
        self.args = args

    def key(self, monomial):
        return tuple(O(lamda(monomial)) for (O, lamda) in self.args)

    def __str__(self):
        from sympy.core import Tuple
        return "ProductOrder" + str(Tuple(*[x[0] for x in self.args]))

    def __eq__(self, other):
        if not isinstance(other, ProductOrder):
            return False
        return self.args == other.args

    def __hash__(self):
        return hash((self.__class__, self.args))

    @property
    def is_global(self):
        if all(o.is_global is True for o, _ in self.args):
            return True
        if all(o.is_global is False for o, _ in self.args):
            return False
        return None

class InverseOrder(MonomialOrder):
    """
    The "inverse" of another monomial order.

    If O is any monomial order, we can construct another monomial order iO
    such that `A >_{iO} B` if and only if `B >_O A`. This is useful for
    constructing local orders.

    Note that many algorithms only work with *global* orders.

    For example, in the inverse lexicographic order on a single variable `x`,
    high powers of `x` count as small:

    >>> from sympy.polys.monomialtools import lex, InverseOrder
    >>> ilex = InverseOrder(lex)
    >>> ilex((5,)) < ilex((0,))
    True
    """

    def __init__(self, O):
        self.O = O

    def __str__(self):
        return "i" + str(self.O)

    def key(self, monomial):
        from sympy.core.compatibility import iterable
        def inv(l):
            if iterable(l):
                return tuple(inv(x) for x in l)
            return -l
        return inv(self.O.key(monomial))

    @property
    def is_global(self):
        if self.O.is_global is True:
            return False
        if self.O.is_global is False:
            return True
        return None

    def __eq__(self, other):
        return isinstance(other, InverseOrder) and other.O == self.O

    def __hash__(self, other):
        return hash((self.__class__, self.O))

lex = LexOrder()
grlex = GradedLexOrder()
grevlex = ReversedGradedLexOrder()
ilex = InverseOrder(lex)
igrlex = InverseOrder(grlex)
igrevlex = InverseOrder(grevlex)

_monomial_key = {
    'lex'      : lex,
    'grlex'    : grlex,
    'grevlex'  : grevlex,
    'ilex'     : ilex,
    'igrlex'   : igrlex,
    'igrevlex' : igrevlex
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
    4. ilex, igrlex, igrevlex - the corresponding inverse orders

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

class _ItemGetter(object):
    """Helper class to return a subsequence of values."""

    def __init__(self, seq):
        self.seq = tuple(seq)

    def __call__(self, m):
        return tuple(m[idx] for idx in self.seq)

    def __eq__(self, other):
        if not isinstance(other, _ItemGetter):
            return False
        return self.seq == other.seq

def build_product_order(arg, gens):
    """
    Build a monomial order on ``gens``.

    ``arg`` should be a tuple of iterables. The first element of each iterable
    should be a string or monomial order (will be passed to monomial_key),
    the others should be subsets of the generators. This function will build
    the corresponding product order.

    For example, build a product of two grlex orders:

    >>> from sympy.polys.monomialtools import grlex, build_product_order
    >>> from sympy.abc import x, y, z, t
    >>> O = build_product_order((("grlex", x, y), ("grlex", z, t)), [x, y, z, t])
    >>> O((1, 2, 3, 4))
    ((3, (1, 2)), (7, (3, 4)))
    """
    gens2idx = {}
    for i, g in enumerate(gens):
        gens2idx[g] = i
    order = []
    for expr in arg:
        name = expr[0]
        var = expr[1:]
        def makelambda(var):
            return _ItemGetter(gens2idx[g] for g in var)
        order.append((monomial_key(name), makelambda(var)))
    return ProductOrder(*order)

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

# TODO cythonize
def monomial_divides(A, B):
    """
    Does there exist a monomial X such that XA == B?

    >>> from sympy.polys.monomialtools import monomial_divides
    >>> monomial_divides((1, 2), (3, 4))
    True
    >>> monomial_divides((1, 2), (0, 2))
    False
    """
    return all(a <= b for a, b in zip(A, B))

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

def monomial_deg(M):
    """
    Returns the total degree of a monomial.

    For example, the total degree of `xy^2` is 3:

    >>> from sympy.polys.monomialtools import monomial_deg
    >>> monomial_deg((1, 2))
    3
    """
    return sum(M)

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
