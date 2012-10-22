from sympy.core import C, Expr, Mul, S, sympify
from sympy.functions.elementary.piecewise import piecewise_fold
from sympy.polys import quo, roots
from sympy.simplify import powsimp

class Product(Expr):
    """Represents unevaluated product.

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
            raise ValueError('Product requires values for lower and upper bounds.')

        obj = Expr.__new__(cls, **assumptions)
        arglist = [sign*function]
        arglist.extend(limits)
        obj._args = tuple(arglist)
        obj.is_commutative = function.is_commutative # limits already checked

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
        from sympy.concrete.summations import _free_symbols

        if self.function.is_zero or self.function == 1:
            return set()
        return _free_symbols(self.function, self.limits)

    @property
    def is_zero(self):
        """A Product is zero only if its term is zero.
        """
        return self.term.is_zero

    @property
    def is_number(self):
        """
        Return True if the Product will result in a number, else False.

        sympy considers anything that will result in a number to have
        is_number == True.

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

    def doit(self, **hints):
        f = g = self.function
        for index, limit in enumerate(self.limits):
            i, a, b = limit
            dif = b - a
            if dif.is_Integer and dif < 0:
                a, b = b, a

            g = self._eval_product(f, (i, a, b))
            if g is None:
                return Product(powsimp(f), *self.limits[index:])
            else:
                f = g

        if hints.get('deep', True):
            return f.doit(**hints)
        else:
            return powsimp(f)

    def _eval_product(self, term, limits):
        from sympy import summation

        (k, a, n) = limits

        if k not in term.free_symbols:
            return term**(n - a + 1)

        if a == n:
            return term.subs(k, a)

        dif = n - a
        if dif.is_Integer:
            return Mul(*[term.subs(k, a + i) for i in xrange(dif  + 1)])

        elif term.is_polynomial(k):
            poly = term.as_poly(k)

            A = B = Q = S.One

            all_roots = roots(poly, multiple=True)

            for r in all_roots:
                A *= C.RisingFactorial(a-r, n-a+1)
                Q *= n - r

            if len(all_roots) < poly.degree():
                arg = quo(poly, Q.as_poly(k))
                B = Product(arg, (k, a, n)).doit()

            return poly.LC()**(n-a+1) * A * B

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
                B = Product(arg, (k, a, n)).doit()
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
                return Product(evaluated, limits)
            else:
                return f

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
    k!
    >>> product(m, (i, 1, k))
    m**k
    >>> product(i, (i, 1, k), (k, 1, n))
    Product(k!, (k, 1, n))

    """

    prod = Product(*args, **kwargs)

    if isinstance(prod, Product):
        return prod.doit(deep=False)
    else:
        return prod
