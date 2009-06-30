
from sympy.core import (Basic, S, C, Add, Mul, Symbol, Equality, Interval,
    sympify, symbols, Wild)
from sympy.functions import factorial
from sympy.solvers import solve

class Sum(Basic):
    """Represents unevaluated summation."""

    def __new__(cls, f, *symbols, **assumptions):
        f = sympify(f)

        if f.is_Number:
            if f is S.NaN:
                return S.NaN
            elif f is S.Zero:
                return S.Zero

        if not symbols:
            limits = f.atoms(Symbol)

            if not limits:
                return f
        else:
            limits = []

            for V in symbols:
                if isinstance(V, Symbol):
                    limits.append(V)
                    continue
                elif isinstance(V, Equality):
                    if isinstance(V.lhs, Symbol):
                        if isinstance(V.rhs, Interval):
                            limits.append((V.lhs, V.rhs.start, V.rhs.end))
                        else:
                            limits.append((V.lhs, V.rhs))

                        continue
                elif isinstance(V, (tuple, list)):
                    if len(V) == 1:
                        if isinstance(V[0], Symbol):
                            limits.append(V[0])
                            continue
                    elif len(V) in (2, 3):
                        if isinstance(V[0], Symbol):
                            limits.append(tuple(map(sympify, V)))
                            continue

                raise ValueError("Invalid summation variable or limits")

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (f, tuple(limits))

        return obj

    @property
    def function(self):
        return self._args[0]

    @property
    def limits(self):
        return self._args[1]

    def doit(self, **hints):
        #if not hints.get('sums', True):
        #    return self
        f = self.function
        for i, a, b in self.limits:
            f = eval_sum(f, (i, a, b))
            if f is None:
                return self
        return f

    def _eval_summation(self, f, x):
        return

    def euler_maclaurin(self, m=0, n=0, eps=0, eval_integral=True):
        """
        Return an Euler-Maclaurin approximation of self, where m is the
        number of leading terms to sum directly and n is the number of
        terms in the tail.

        With m = n = 0, this is simply the corresponding integral
        plus a first-order endpoint correction.

        Returns (s, e) where s is the Euler-Maclaurin approximation
        and e is the estimated error (taken to be the magnitude of
        the first omitted term in the tail):

            >>> k = Symbol('k')
            >>> Sum(1/k, (k, 2, 5)).doit().evalf()
            1.28333333333333
            >>> s, e = Sum(1/k, (k, 2, 5)).euler_maclaurin()
            >>> s
            7/20 - log(2) + log(5)
            >>> from sympy import sstr
            >>> print sstr((s.evalf(), e.evalf()), full_prec=True)
            (1.26629073187416, 0.0175000000000000)

        The endpoints may be symbolic:

            >>> k, a, b = symbols('kab')
            >>> s, e = Sum(1/k, (k, a, b)).euler_maclaurin()
            >>> s
            -log(a) + log(b) + 1/(2*a) + 1/(2*b)
            >>> e
            abs(-1/(12*b**2) + 1/(12*a**2))

        If the function is a polynomial of degree at most 2n+1, the
        Euler-Maclaurin formula becomes exact (and e = 0 is returned):

            >>> Sum(k, (k, 2, b)).euler_maclaurin()
            (-1 + b/2 + b**2/2, 0)
            >>> Sum(k, (k, 2, b)).doit()
            -1 + b/2 + b**2/2

        With a nonzero eps specified, the summation is ended
        as soon as the remainder term is less than the epsilon.
        """
        m = int(m)
        n = int(n)
        f = self.function
        assert len(self.limits) == 1
        i, a, b = self.limits[0]
        s = S.Zero
        if m:
            for k in range(m):
                term = f.subs(i, a+k)
                if (eps and term and abs(term.evalf(3)) < eps):
                    return s, abs(term)
                s += term
            a += m
        x = Symbol('x', dummy=True)
        I = C.Integral(f.subs(i, x), (x, a, b))
        if eval_integral:
            I = I.doit()
        s += I
        def fpoint(expr):
            if b is S.Infinity:
                return expr.subs(i, a), 0
            return expr.subs(i, a), expr.subs(i, b)
        fa, fb = fpoint(f)
        iterm = (fa + fb)/2
        g = f.diff(i)
        for k in xrange(1, n+2):
            ga, gb = fpoint(g)
            term = C.bernoulli(2*k)/C.Factorial(2*k)*(gb-ga)
            if (eps and term and abs(term.evalf(3)) < eps) or (k > n):
                break
            s += term
            g = g.diff(i, 2)
        return s + iterm, abs(term)

    def _eval_subs(self, old, new):
        newargs = (self.args[1][0][0], self.args[1][0][1].subs(old,new),
                   self.args[1][0][2].subs(old,new))
        return Sum(self.args[0].subs(old, new), newargs)


def sum(*args, **kwargs):
    summation = Sum(*args, **kwargs)

    if isinstance(summation, Sum):
        return summation.doit()
    else:
        return summation


def getab(expr):
    cls = expr.__class__
    return cls(expr.args[0]), cls(*expr.args[1:])

def telescopic_direct(L, R, n, (i, a, b)):
    '''Returns the direct summation of the terms of a telescopic sum

    L is the term with lower index
    R is the term with higher index
    n difference between the indexes of L and R

    For example:

    >>> k,a,b = symbols('kab')
    >>> telescopic_direct(1/k, -1/(k+2), 2, (k, a, b))
    1/a + 1/(1 + a) - 1/(1 + b) - 1/(2 + b)
    '''
    s = 0
    for m in xrange(n):
        s += L.subs(i,a+m) + R.subs(i,b-m)
    return s

def telescopic(L, R, (i, a, b)):
    '''Tries to perform the summation using the telescopic property

    return None if not possible
    '''
    if L.is_Add or R.is_Add:
        return None
    s = None
    #First we try to solve using match
    #Maybe this should go inside solve
    k = Wild("k")
    sol = (-R).match(L.subs(i, i + k))
    if sol and k in sol:
        if L.subs(i,i + sol[k]) == -R:
            #sometimes match fail(f(x+2).match(-f(x+k))->{k: -2 - 2x}))
            s = sol[k]
    #Then we try to solve using solve
    if not s or not s.is_Integer:
        m = Symbol("m")
        try:
            s = solve(L.subs(i, i + m) + R, m)[0]
        except ValueError:
            pass
    if s and s.is_Integer:
        if s < 0:
            return telescopic_direct(R, L, abs(s), (i, a, b))
        elif s > 0:
            return telescopic_direct(L, R, s, (i, a, b))
    return None

def eval_sum(f, (i, a, b)):
    if not f.has(i):
        return f*(b-a+1)
    definite = a.is_Integer and b.is_Integer
    # Doing it directly may be faster if there are very few terms.
    if definite and (b-a < 100):
        return eval_sum_direct(f, (i, a, b))
    # Try to do it symbolically. Even when the number of terms is known,
    # this can save time when b-a is big.
    # We should try to transform to partial fractions
    value = eval_sum_symbolic(f.expand(), (i, a, b))
    if value is not None:
        return value
    # Do it directly
    if definite:
        return eval_sum_direct(f, (i, a, b))

def eval_sum_symbolic(f, (i, a, b)):
    if not f.has(i):
        return f*(b-a+1)
    # Linearity
    if f.is_Mul:
        L, R = getab(f)
        if not L.has(i):
            sR = eval_sum_symbolic(R, (i, a, b))
            if sR: return L*sR
        if not R.has(i):
            sL = eval_sum_symbolic(L, (i, a, b))
            if sL: return R*sL
    if f.is_Add:
        L, R = getab(f)
        lrsum = telescopic(L, R, (i, a, b))
        if lrsum: return lrsum
        lsum = eval_sum_symbolic(L, (i, a, b))
        rsum = eval_sum_symbolic(R, (i, a, b))
        if None not in (lsum, rsum):
            return lsum + rsum
    # Polynomial terms with Faulhaber's formula
    p = C.Wild('p')
    e = f.match(i**p)
    if e != None:
        c = p.subs(e)
        B = C.bernoulli
        if c.is_integer and c >= 0:
            s = (B(c+1, b+1) - B(c+1, a))/(c+1)
            return s.expand()
    # Geometric terms
    c1 = C.Wild('c1', exclude=[i])
    c2 = C.Wild('c2', exclude=[i])
    c3 = C.Wild('c3', exclude=[i])
    e = f.match(c1**(c2*i+c3))
    if e is not None:
        c1 = c1.subs(e)
        c2 = c2.subs(e)
        c3 = c3.subs(e)
        # TODO: more general limit handling
        return c1**c3 * (c1**(a*c2) - c1**(c2+b*c2)) / (1 - c1**c2)
    return None

def eval_sum_direct(expr, (i, a, b)):
    s = S.Zero
    if expr.has(i):
        for j in xrange(a, b+1):
            s += expr.subs(i, j)
    else:
        for j in xrange(a, b+1):
            s += expr
    return s
