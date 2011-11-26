from sympy.core import Expr, S, C, Mul, sympify, Symbol, Tuple
from sympy.core.compatibility import is_sequence
from sympy.polys import quo, roots
from sympy.simplify import powsimp

class Product(Expr):
    """Represents unevaluated product.

    """

    def __new__(cls, term, *symbols, **assumptions):
        term = sympify(term)

        if term is S.NaN:
            return S.NaN

        if len(symbols) == 1:
            symbol = symbols[0]

            if isinstance(symbol, C.Equality):
                k = symbol.lhs
                a = symbol.rhs.start
                n = symbol.rhs.end
            elif is_sequence(symbol):
                k, a, n = symbol
            else:
                raise ValueError("Invalid arguments")

            k, a, n = map(sympify, (k, a, n))

        else:
            raise NotImplementedError

        obj = Expr.__new__(cls, **assumptions)
        obj._args = (term, Tuple(k, a, n))

        return obj

    @property
    def term(self):
        return self._args[0]
    function = term

    @property
    def index(self):
        return self._args[1][0]

    @property
    def lower(self):
        return self._args[1][1]

    @property
    def upper(self):
        return self._args[1][2]

    @property
    def limits(self):
        return (self._args[1],)

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
        term = self.term
        if term == 0:
            return S.Zero
        elif term == 1:
            return S.One
        lower = self.lower
        upper = self.upper
        if hints.get('deep', True):
            term = term.doit(**hints)
            lower = lower.doit(**hints)
            upper = upper.doit(**hints)
        dif = upper - lower
        if dif.is_Number and dif < 0:
            upper, lower = lower, upper

        prod = self._eval_product(lower, upper, term)

        if prod is not None:
            return powsimp(prod)
        else:
            return self

    def _eval_product(self, a, n, term):
        from sympy import summation

        k = self.index

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

            p = self._eval_product(a, n, p)
            q = self._eval_product(a, n, q)

            return p / q

        elif term.is_Mul:
            exclude, include = [], []

            for t in term.args:
                p = self._eval_product(a, n, t)

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
                p = self._eval_product(a, n, term.base)

                if p is not None:
                    return p**term.exp

def product(*args, **kwargs):
    prod = Product(*args, **kwargs)

    if isinstance(prod, Product):
        return prod.doit(deep=False)
    else:
        return prod

