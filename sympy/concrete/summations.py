
from sympy.core import Basic, S, C, Add, Mul, Symbol, Equality, Interval, sympify
from sympy.core.methods import NoRelMeths, ArithMeths

class Sum(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated summation."""

    precedence = Basic.Apply_precedence

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

    def tostr(self, level=0):
        L = ', '.join([ str(L) for L in self.limits ])
        return 'Sum(%s, %s)' % (self.function.tostr(), L)

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

def sum(*args, **kwargs):
    summation = Sum(*args, **kwargs)

    if isinstance(summation, Sum):
        return summation.doit()
    else:
        return summation


def getab(expr):
    cls = expr.__class__
    return cls(expr.args[0]), cls(*expr.args[1:])

def eval_sum(f, (i, a, b)):
    if not f.has(i):
        return f*(b-a+1)
    definite = a.is_Integer and b.is_Integer
    # Doing it directly may be faster if there are very few terms.
    if definite and (b-a < 100):
        return eval_sum_direct(f, (i, a, b))
    # Try to do it symbolically. Even when the number of terms is known,
    # this can save time when b-a is big.
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
        if not L.has(i): return L*eval_sum_symbolic(R, (i, a, b))
        if not R.has(i): return R*eval_sum_symbolic(L, (i, a, b))
    if f.is_Add:
        L, R = getab(f)
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
