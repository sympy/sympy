from sympy.core import Basic, C, Rational, Add, Mul, Pow, Symbol, Wild, oo
from sympy.core import sympify
#from sympy.specfun import rising_factorial, factorial, factorial_simplify
#from sympy.specfun.factorials import unfac
#from sympy.specfun import bernoulli
#from sympy.simplify import powsimp

def ispoly(expr, var):
    return False

def getab(expr):
    cls = expr.__class__
    return cls(expr.args[0]), cls(*expr.args[1:])

def indexsymbol(a):
    if isinstance(a, Symbol):
        return Symbol(a.name, integer=True)
    else:
        return sympify(a)

class _BigOperator(Basic):

    def __new__(cls, f, (i, a, b)):
        self = Basic.__new__(cls)
        assert isinstance(i, Symbol)
        self.i = i
        self.f = sympify(f)
        self.a = indexsymbol(a)
        self.b = indexsymbol(b)
        return self.eval()


class Sum2(_BigOperator):
    """
    Symbolic summation with a variable number of terms

    Sum2(f, (i, a, b)) represents \sum_{i=a}^b f(i)
    """

    def reindex(self, a):
        """Re-index the sum to start at a new lower index a."""
        diff = self.a - a
        b = self.b - diff
        f = self.f.subs(self.i, self.i + diff)
        return Sum2(f, (self.i, a, b))

    def split(self, n):
        """Split into two sums, the first with n terms."""
        f, i, a, b = self.f, self.i, self.a, self.b
        return Sum2(f, (i, a, a+n-1)) + Sum2(f, (i, a+n, b))

    def eval(self):
        f, i, a, b = self.f, self.i, self.a, self.b

        # Exploit the linearity of the sum
        if not f.has(i):
            return f*(b-a+1)
        if f.is_Mul:
            L, R = getab(f)
            if not L.has(i): return L*Sum2(R, (i, a, b))
            if not R.has(i): return R*Sum2(L, (i, a, b))
        if f.is_Add:
            L, R = getab(f)
            lsum = Sum2(L, (i,a,b))
            rsum = Sum2(R, (i,a,b))
            if not isinstance(lsum, Sum2) and not isinstance(rsum, Sum2):
                return lsum + rsum

        # Polynomial terms with Faulhaber's formula
        if f == i:
            f = Pow(i, 1, evaluate=False) # TODO: match should handle this
        p = Wild('p')
        e = f.match(i**p)
        if e != None:
            c = p.subs(e)
            B = C.bernoulli
            if c.is_integer and c >= 0:
                s = (B(c+1, b+1) - B(c+1, a))/(c+1)
                return s.expand()

        # Geometric terms
        if f.is_Pow:
            r, k = f.args
            if not r.has(i) and k == i:
                # TODO: Pow should be able to simplify x**oo depending
                # on whether |x| < 1 or |x| > 1 for non-rational x
                if b == oo and isinstance(r, Rational) and abs(r) < 1:
                    return r**a / (1-r)
                else:
                    return (r**a - r**(b+1)) / (1-r)

        # Should nothing else works, use brute force if possible
        if a.is_Integer and b.is_Integer:
            s = 0
            for j in range(a, b+1):
                s += f.subs(i, j)
            return s

        return self

    def subs(self, x, y):
        if x == self.b:
            return Sum2(self.f, (self.i, self.a, y))
        return self

'''
class Product(_BigOperator):
    """
    Symbolic product with a variable number of factors

    Product(f, (i, a, b)) represents \prod_{i=a}^b f(i)
    """

    def __repr__(self):
        return "Product(%r, (%r, %r, %r))" % (self.f, self.i, self.a, self.b)

    __str__ = __repr__

    def eval(self):
        # Simplify sub-products
        p = self._eval()
        if isinstance(p, Product):
            return self
        else:
            return powsimp(factorial_simplify(p))
        return p

    def _eval(self):
        f, i, a, b = self.f, self.i, self.a, self.b

        if not f.has(i):
            return f**(b-a+1)

        if f.is_Mul:
            L, R = getab(f)
            lp = Product(L, (i, a, b))
            rp = Product(R, (i, a, b))
            if not (isinstance(lp, Product) and isinstance(rp, Product)):
                return lp * rp

        if f.is_Pow:
            base, exp = f.args
            if not base.has(i):
                s = Sum(exp, (i, a, b))
                if not isinstance(s, Sum):
                    return base ** s
            elif not exp.has(i):
                p = Product(base, (i, a, b))
                if not isinstance(p, Product):
                    return p ** exp

        # Linear functions
        if f == i:
            return rising_factorial(a, b-a+1)
        #if ispoly(f, i):
        p = Wild('p')
        q = Wild('q')
        e = f.match(p+q*i)
        if e != None:
            pp = p.subs(e)
            qq = q.subs(e)
            if not pp.has(i) and not qq.has(i):
                r = qq**(b-a+1) * unfac(b+pp/qq) / unfac(a+pp/qq-1)
                return r

        # Given a more complicated rational expression, try to factor
        # it into linear functions
        if f.is_Add:
            try:
                num, den = fraction(together(f))
                g = factor(num) / factor(den)
                p = Product(g, (i, a, b))
                if not isinstance(p, Product):
                    return p
            except PolynomialException:
                pass

        # Brute force
        if a.is_Integer and b.is_Integer:
            p = 1
            for j in range(a, b+1):
                p *= f.subs(i, j)
            return p

        return self
'''
