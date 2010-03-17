from sympy.core import Basic, S, C, Mul, sympify

from sympy.polys import quo, roots
from sympy.simplify import powsimp

class Product(Basic):
    """Represents unevaluated product.

    """

    def __new__(cls, term, *symbols, **assumptions):
        term = sympify(term)

        if term.is_Number:
            if term is S.NaN:
                return S.NaN
            elif term is S.Infinity:
                return S.NaN
            elif term is S.NegativeInfinity:
                return S.NaN
            elif term is S.Zero:
                return S.Zero
            elif term is S.One:
                return S.One

        if len(symbols) == 1:
            symbol = symbols[0]

            if isinstance(symbol, C.Equality):
                k = symbol.lhs
                a = symbol.rhs.start
                n = symbol.rhs.end
            elif isinstance(symbol, (tuple, list)):
                k, a, n = symbol
            else:
                raise ValueError("Invalid arguments")

            k, a, n = map(sympify, (k, a, n))

            if isinstance(a, C.Number) and isinstance(n, C.Number):
                return Mul(*[term.subs(k, i) for i in xrange(int(a), int(n)+1)])
        else:
            raise NotImplementedError

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (term, k, a, n)

        return obj

    @property
    def term(self):
        return self._args[0]

    @property
    def index(self):
        return self._args[1]

    @property
    def lower(self):
        return self._args[2]

    @property
    def upper(self):
        return self._args[3]

    def doit(self, **hints):
        term = self.term
        lower = self.lower
        upper = self.upper
        if hints.get('deep', True):
            term = term.doit(**hints)
            lower = lower.doit(**hints)
            upper = upper.doit(**hints)

        prod = self._eval_product(lower, upper, term)

        if prod is not None:
            return powsimp(prod)
        else:
            return self

    def _eval_product(self, a, n, term):
        from sympy import sum, Sum
        k = self.index

        if not term.has(k):
            return term**(n-a+1)
        elif term.is_polynomial(k):
            poly = term.as_poly(k)

            A = B = Q = S.One
            C_= poly.LC()

            all_roots = roots(poly, multiple=True)

            for r in all_roots:
                A *= C.RisingFactorial(a-r, n-a+1)
                Q *= n - r

            if len(all_roots) < poly.degree():
                B = Product(quo(poly, Q.as_poly(k)), (k, a, n))

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
                A, B = Mul(*exclude), Mul(*include)
                return A * Product(B, (k, a, n))
        elif term.is_Pow:
            if not term.base.has(k):
                s = sum(term.exp, (k, a, n))

                if not isinstance(s, Sum):
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

