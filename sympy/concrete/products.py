
from sympy.core import Basic, S, C, Add, Mul, Symbol, sympify
from sympy.core.methods import NoRelMeths, ArithMeths

from sympy.polynomials import quo, roots
from sympy.simplify import powsimp

class Product(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated product.

    """

    precedence = Basic.Apply_precedence

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

    def tostr(self, level=0):
        r = 'Product(%s)' % ', '.join([a.tostr() for a in self._args])

        if self.precedence <= level:
            return '(%s)' % (r)
        else:
            return r

    def doit(self):
        prod = self._eval_product()

        if prod is not None:
            return powsimp(prod)
        else:
            return self

    def _eval_product(self, term=None):
        k = self.index
        a = self.lower
        n = self.upper

        if term is None:
            term = self.term

        if not term.has(k):
            return term**(n-a+1)
        elif term.is_polynomial(k):
            poly = term.as_polynomial(k)

            A = B = Q = S.One
            C_= poly.leading_coeff()

            all_roots = roots(poly)

            for r in all_roots:
                A *= C.RisingFactorial(a-r, n-a+1)
                Q *= n - r

            if len(all_roots) < poly.degree():
                B = Product(quo(poly, Q, k), (k, a, n))

            return C_**(n-a+1) * A * B
        elif term.is_Add:
            p, q = term.as_numer_denom()

            p = self._eval_product(p)
            q = self._eval_product(q)

            return p / q
        elif term.is_Mul:
            exclude, include = [], []

            for t in term.args:
                p = self._eval_product(t)

                if p is not None:
                    exclude.append(p)
                else:
                    include.append(p)

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
                p = self._eval_product(term.base)

                if p is not None:
                    return p**term.exp

def product(*args, **kwargs):
    prod = Product(*args, **kwargs)

    if isinstance(prod, Product):
        return prod.doit()
    else:
        return prod
