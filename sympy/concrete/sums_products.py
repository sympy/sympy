"""
from sympy import *
from sympy.specfun import rising_factorial, factorial, factorial_simplify
from sympy.specfun.factorials import unfac
from sympy.specfun import bernoulli
from sympy.polynomials import factor, PolynomialException
from sympy.simplify import powsimp

def ispoly(expr, var):
    return False

def getab(expr):
    cls = expr.__class__
    return cls(expr[0]), cls(*expr[1:])

def indexsymbol(a):
    if isinstance(a, Symbol):
        return Symbol(a.name, integer=True)
    else:
        return Basic.sympify(a)

class _BigOperator(Basic):

    def __new__(cls, f, (i, a, b)):
        self = Basic.__new__(cls)
        assert isinstance(i, Symbol)
        self.i = i
        self.f = self.sympify(f)
        self.a = indexsymbol(a)
        self.b = indexsymbol(b)
        return self.eval()


class Sum(_BigOperator):

    Symbolic summation with a variable number of terms

    Sum(f, (i, a, b)) represents \sum_{i=a}^b f(i)


    def __repr__(self):
        return "Sum(%r, (%r, %r, %r))" % (self.f, self.i, self.a, self.b)

    __str__ = __repr__

    def reindex(self, a):
        Reindex the sum to start at a new lower index a.
        diff = self.a - a
        b = self.b - diff
        f = self.f.subs(self.i, self.i + diff)
        return Sum(f, (self.i, a, b))

    def split(self, n):
        Split into two sums, the first with n terms.
        f, i, a, b = self.f, self.i, self.a, self.b
        return Sum(f, (i, a, a+n-1)) + Sum(f, (i, a+n, b))

    def eval(self):
        f, i, a, b = self.f, self.i, self.a, self.b

        # Exploit the linearity of the sum
        if not f.has(i):
            return f*(b-a+1)
        if isinstance(f, Mul):
            L, R = getab(f)
            if not L.has(i): return L*Sum(R, (i, a, b))
            if not R.has(i): return R*Sum(L, (i, a, b))
        if isinstance(f, Add):
            L, R = getab(f)
            lsum = Sum(L, (i,a,b))
            rsum = Sum(R, (i,a,b))
            if not isinstance(lsum, Sum) and not isinstance(rsum, Sum):
                return lsum + rsum

        # Polynomial terms with Faulhaber's formula
        if f == i:
            f = Pow(i, 1, evaluate=False) # TODO: match should handle this
        p = Wild('p')
        e = f.match(i**p)
        if e != None:
            c = p.subs_dict(e)
            if c.is_integer and c >= 0:
                s = (bernoulli(c+1, b+1)-bernoulli(c+1, a))/(c+1)
                return s.expand()

        # Geometric terms
        if isinstance(f, Pow):
            r, k = f[:]
            if not r.has(i) and k == i:
                # TODO: Pow should be able to simplify x**oo depending
                # on whether |x| < 1 or |x| > 1 for non-rational x
                if b == oo and isinstance(r, Rational) and abs(r) < 1:
                    return r**a / (1-r)
                else:
                    return (r**a - r**(b+1)) / (1-r)

        # Should nothing else works, use brute force if possible
        if isinstance(a, Rational) and a.is_integer and \
           isinstance(b, Rational) and b.is_integer:
            s = 0
            for j in range(a, b+1):
                s += f.subs(i, j)
            return s

        return self

    def euler_maclaurin(self, n=0):

        Return n-th order Euler-Maclaurin approximation of self.

        The 0-th order approximation is simply the corresponding
        integral

        f, i, a, b = self.f, self.i, self.a, self.b
        x = Symbol('x')
        s = integrate(f.subs(i, x), x==[a,b])
        if n > 0:
            s += (f.subs(i, a) + f.subs(i, b))/2
        for k in range(1, n):
            g = f.diff(i, 2*k-1)
            s += bernoulli(2*k)/factorial(2*k)*(g.subs(i,b)-g.subs(i,a))
        return s


class Product(_BigOperator):

    Symbolic product with a variable number of factors

    Product(f, (i, a, b)) represents \prod_{i=a}^b f(i)


    def __repr__(self):
        return "Product(%r, (%r, %r, %r))" % (self.f, self.i, self.a, self.b)

    __str__ = __repr__

    def eval(self):
        # Simplify subproducts
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

        if isinstance(f, Mul):
            L, R = getab(f)
            lp = Product(L, (i, a, b))
            rp = Product(R, (i, a, b))
            if not (isinstance(lp, Product) and isinstance(rp, Product)):
                return lp * rp

        if isinstance(f, Pow):
            base, exp = f[:]
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
            pp = p.subs_dict(e)
            qq = q.subs_dict(e)
            if not pp.has(i) and not qq.has(i):
                r = qq**(b-a+1) * unfac(b+pp/qq) / unfac(a+pp/qq-1)
                return r

        # Given a more complicated rational expression, try to factor
        # it into linear functions
        if isinstance(f, Add):
            try:
                num, den = fraction(together(f))
                g = factor(num) / factor(den)
                p = Product(g, (i, a, b))
                if not isinstance(p, Product):
                    return p
            except PolynomialException:
                pass

        # Brute force
        if isinstance(a, Rational) and a.is_integer and \
           isinstance(b, Rational) and b.is_integer:
            p = 1
            for j in range(a, b+1):
                p *= f.subs(i, j)
            return p

        return self
"""