from sympy.core.containers import Tuple
from sympy.core.core import C
from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.functions.elementary.piecewise import piecewise_fold
from sympy.polys import quo, roots
from sympy.simplify import powsimp


class Product(Expr):
    r"""Represents unevaluated products.

    ``Product`` represents a finite or infinite product, with the first
    argument being the general form of terms in the series, and the second
    argument being ``(dummy_variable, start, end)``, with ``dummy_variable``
    taking all integer values from ``start`` through ``end``. In accordance
    with long-standing mathematical convention, the end term is included in
    the product.

    Finite products
    ===============

    For finite products (and products with symbolic limits assumed to be finite)
    we follow the analogue of the summation convention described by Karr [1],
    especially definition 3 of section 1.4. The product:

    .. math::

        \prod_{m \leq i < n} f(i)

    has *the obvious meaning* for `m < n`, namely:

    .. math::

        \prod_{m \leq i < n} f(i) = f(m) f(m+1) \cdot \ldots \cdot f(n-2) f(n-1)

    with the upper limit value `f(n)` excluded. The product over an empty set is
    one if and only if `m = n`:

    .. math::

        \prod_{m \leq i < n} f(i) = 1  \quad \mathrm{for} \quad  m = n

    Finally, for all other products over empty sets we assume the following
    definition:

    .. math::

        \prod_{m \leq i < n} f(i) = \frac{1}{\prod_{n \leq i < m} f(i)}  \quad \mathrm{for} \quad  m > n

    It is important to note that above we define all products with the upper
    limit being exclusive. This is in contrast to the usual mathematical notation,
    but does not affect the product convention. Indeed we have:

    .. math::

        \prod_{m \leq i < n} f(i) = \prod_{i = m}^{n - 1} f(i)

    where the difference in notation is intentional to emphasize the meaning,
    with limits typeset on the top being inclusive.

    Examples
    ========

    >>> from sympy.abc import a, b, i, k, m, n, x
    >>> from sympy import Product, factorial, oo
    >>> Product(k,(k,1,m))
    Product(k, (k, 1, m))
    >>> Product(k,(k,1,m)).doit()
    factorial(m)
    >>> Product(k**2,(k,1,m))
    Product(k**2, (k, 1, m))
    >>> Product(k**2,(k,1,m)).doit()
    (factorial(m))**2

    Wallis' product for pi:

    >>> W = Product(2*i/(2*i-1) * 2*i/(2*i+1), (i, 1, oo))
    >>> W
    Product(4*i**2/((2*i - 1)*(2*i + 1)), (i, 1, oo))

    Direct computation currently fails:

    >>> W.doit()
    nan

    But we can approach the infinite product by a limit of finite products:

    >>> from sympy import limit
    >>> W2 = Product(2*i/(2*i-1)*2*i/(2*i+1), (i, 1, n))
    >>> W2
    Product(4*i**2/((2*i - 1)*(2*i + 1)), (i, 1, n))
    >>> W2e = W2.doit()
    >>> W2e
    2**(-2*n)*4**n*(factorial(n))**2/(RisingFactorial(1/2, n)*RisingFactorial(3/2, n))
    >>> limit(W2e, n, oo)
    pi/2

    By the same formula we can compute sin(pi/2):

    >>> from sympy import pi, gamma, simplify
    >>> P = pi * x * Product(1 - x**2/k**2,(k,1,n))
    >>> P = P.subs(x, pi/2)
    >>> P
    pi**2*Product(1 - pi**2/(4*k**2), (k, 1, n))/2
    >>> Pe = P.doit()
    >>> Pe
    pi**2*RisingFactorial(1 + pi/2, n)*RisingFactorial(-pi/2 + 1, n)/(2*(factorial(n))**2)
    >>> Pe = Pe.rewrite(gamma)
    >>> Pe
    pi**2*gamma(n + 1 + pi/2)*gamma(n - pi/2 + 1)/(2*gamma(1 + pi/2)*gamma(-pi/2 + 1)*gamma(n + 1)**2)
    >>> Pe = simplify(Pe)
    >>> Pe
    sin(pi**2/2)*gamma(n + 1 + pi/2)*gamma(n - pi/2 + 1)/gamma(n + 1)**2
    >>> limit(Pe, n, oo)
    sin(pi**2/2)

    Products with the lower limit being larger than the upper one:

    >>> Product(1/i, (i, 6, 1)).doit()
    120
    >>> Product(i, (i, 2, 5)).doit()
    120

    The empty product:

    >>> Product(i, (i, n, n-1)).doit()
    1

    An example showing that the symbolic result of a product is still
    valid for seemingly nonsensical values of the limits. Then the Karr
    convention allows us to give a perfectly valid interpretation to
    those products by interchanging the limits according to the above rules:

    >>> P = Product(2, (i, 10, n)).doit()
    >>> P
    2**(n - 9)
    >>> P.subs(n, 5)
    1/16
    >>> Product(2, (i, 10, 5)).doit()
    1/16
    >>> 1/Product(2, (i, 6, 9)).doit()
    1/16

    An explicit example of the Karr summation convention applied to products:

    >>> P1 = Product(x, (i, a, b)).doit()
    >>> P1
    x**(-a + b + 1)
    >>> P2 = Product(x, (i, b+1, a-1)).doit()
    >>> P2
    x**(a - b - 1)
    >>> simplify(P1 * P2)
    1

    And another one:

    >>> P1 = Product(i, (i, b, a)).doit()
    >>> P1
    RisingFactorial(b, a - b + 1)
    >>> P2 = Product(i, (i, a+1, b-1)).doit()
    >>> P2
    RisingFactorial(a + 1, -a + b - 1)
    >>> P1 * P2
    RisingFactorial(b, a - b + 1)*RisingFactorial(a + 1, -a + b - 1)
    >>> simplify(P1 * P2)
    1

    See Also
    ========

    Sum, summation
    product

    References
    ==========

    .. [1] Michael Karr, "Summation in Finite Terms", Journal of the ACM,
           Volume 28 Issue 2, April 1981, Pages 305-350
           http://dl.acm.org/citation.cfm?doid=322248.322255
    .. [2] http://en.wikipedia.org/wiki/Multiplication#Capital_Pi_notation
    .. [3] http://en.wikipedia.org/wiki/Empty_product
    """

    __slots__ = ['is_commutative']

    def __new__(cls, function, *symbols, **assumptions):
        from sympy.integrals.integrals import _process_limits

        # Any embedded piecewise functions need to be brought out to the
        # top level so that integration can go into piecewise mode at the
        # earliest possible moment.
        function = piecewise_fold(sympify(function))

        if function is S.NaN:
            return S.NaN

        if not symbols:
            raise ValueError("Product variables must be given")

        limits, sign = _process_limits(*symbols)

        # Only limits with lower and upper bounds are supported; the indefinite
        # Product is not supported
        if any(len(l) != 3 or None in l for l in limits):
            raise ValueError(
                'Product requires values for lower and upper bounds.')

        obj = Expr.__new__(cls, **assumptions)
        arglist = [sign*function]
        arglist.extend(limits)
        obj._args = tuple(arglist)
        obj.is_commutative = function.is_commutative  # limits already checked

        return obj

    @property
    def term(self):
        return self._args[0]
    function = term

    @property
    def limits(self):
        return self._args[1:]

    @property
    def variables(self):
        """Return a list of the product variables

        >>> from sympy import Product
        >>> from sympy.abc import x, i
        >>> Product(x**i, (i, 1, 3)).variables
        [i]
        """
        return [l[0] for l in self.limits]

    @property
    def free_symbols(self):
        """
        This method returns the symbols that will affect the value of
        the Product when evaluated. This is useful if one is trying to
        determine whether a product depends on a certain symbol or not.

        >>> from sympy import Product
        >>> from sympy.abc import x, y
        >>> Product(x, (x, y, 1)).free_symbols
        set([y])
        """
        from sympy.integrals.integrals import _free_symbols
        if self.function.is_zero or self.function == 1:
            return set()
        return _free_symbols(self)

    @property
    def is_zero(self):
        """A Product is zero only if its term is zero.
        """
        return self.term.is_zero

    @property
    def is_number(self):
        """
        Return True if the Product will result in a number, else False.

        Examples
        ========

        >>> from sympy import log, Product
        >>> from sympy.abc import x, y, z
        >>> log(2).is_number
        True
        >>> Product(x, (x, 1, 2)).is_number
        True
        >>> Product(y, (x, 1, 2)).is_number
        False
        >>> Product(1, (x, y, z)).is_number
        True
        >>> Product(2, (x, y, z)).is_number
        False
        """

        return self.function.is_zero or self.function == 1 or not self.free_symbols

    def as_dummy(self):
        from sympy.integrals.integrals import _as_dummy
        return _as_dummy(self)

    def doit(self, **hints):
        f = self.function

        for index, limit in enumerate(self.limits):
            i, a, b = limit
            dif = b - a
            if dif.is_Integer and dif < 0:
                a, b = b + 1, a - 1
                f = 1 / f

            g = self._eval_product(f, (i, a, b))
            if g is None:
                return self.func(powsimp(f), *self.limits[index:])
            else:
                f = g

        if hints.get('deep', True):
            return f.doit(**hints)
        else:
            return powsimp(f)

    def _eval_adjoint(self):
        if self.is_commutative:
            return self.func(self.function.adjoint(), *self.limits)
        return None

    def _eval_conjugate(self):
        return self.func(self.function.conjugate(), *self.limits)

    def _eval_product(self, term, limits):
        from sympy.concrete.delta import deltaproduct, _has_simple_delta
        from sympy.concrete.summations import summation
        from sympy.functions import KroneckerDelta

        (k, a, n) = limits

        if k not in term.free_symbols:
            return term**(n - a + 1)

        if a == n:
            return term.subs(k, a)

        if term.has(KroneckerDelta) and _has_simple_delta(term, limits[0]):
            return deltaproduct(term, limits)

        dif = n - a
        if dif.is_Integer:
            return Mul(*[term.subs(k, a + i) for i in range(dif + 1)])

        elif term.is_polynomial(k):
            poly = term.as_poly(k)

            A = B = Q = S.One

            all_roots = roots(poly, multiple=True)

            for r in all_roots:
                A *= C.RisingFactorial(a - r, n - a + 1)
                Q *= n - r

            if len(all_roots) < poly.degree():
                arg = quo(poly, Q.as_poly(k))
                B = self.func(arg, (k, a, n)).doit()

            return poly.LC()**(n - a + 1) * A * B

        elif term.is_Add:
            p, q = term.as_numer_denom()

            p = self._eval_product(p, (k, a, n))
            q = self._eval_product(q, (k, a, n))

            return p / q

        elif term.is_Mul:
            exclude, include = [], []

            for t in term.args:
                p = self._eval_product(t, (k, a, n))

                if p is not None:
                    exclude.append(p)
                else:
                    include.append(t)

            if not exclude:
                return None
            else:
                arg = term._new_rawargs(*include)
                A = Mul(*exclude)
                B = self.func(arg, (k, a, n)).doit()
                return A * B

        elif term.is_Pow:
            if not term.base.has(k):
                s = summation(term.exp, (k, a, n))

                return term.base**s
            elif not term.exp.has(k):
                p = self._eval_product(term.base, (k, a, n))

                if p is not None:
                    return p**term.exp

        elif isinstance(term, Product):
            evaluated = term.doit()
            f = self._eval_product(evaluated, limits)
            if f is None:
                return self.func(evaluated, limits)
            else:
                return f

    def _eval_transpose(self):
        if self.is_commutative:
            return self.func(self.function.transpose(), *self.limits)
        return None


    def _eval_subs(self, old, new):
        from sympy.integrals.integrals import _eval_subs
        return _eval_subs(self, old, new)


def product(*args, **kwargs):
    r"""
    Compute the product.

    The notation for symbols is similiar to the notation used in Sum or
    Integral. product(f, (i, a, b)) computes the product of f with
    respect to i from a to b, i.e.,

    ::

                                     b
                                   _____
        product(f(n), (i, a, b)) = |   | f(n)
                                   |   |
                                   i = a

    If it cannot compute the product, it returns an unevaluated Product object.
    Repeated products can be computed by introducing additional symbols tuples::

    >>> from sympy import product, symbols
    >>> i, n, m, k = symbols('i n m k', integer=True)

    >>> product(i, (i, 1, k))
    factorial(k)
    >>> product(m, (i, 1, k))
    m**k
    >>> product(i, (i, 1, k), (k, 1, n))
    Product(factorial(k), (k, 1, n))

    """

    prod = Product(*args, **kwargs)

    if isinstance(prod, Product):
        return prod.doit(deep=False)
    else:
        return prod
