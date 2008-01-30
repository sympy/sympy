
from sympy.core import Basic, S, C, Add, Mul, Symbol, Equality, Interval, sympify
from sympy.core.methods import NoRelMeths, ArithMeths

class Sum(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated summation."""

    precedence = Basic.Apply_precedence

    def __new__(cls, f, *symbols, **assumptions):
        f = Basic.sympify(f)

        if isinstance(f, C.Number):
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
    definite = isinstance(a, C.Integer) and isinstance(b, C.Integer)
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
    if isinstance(f, C.Mul):
        L, R = getab(f)
        if not L.has(i): return L*eval_sum_symbolic(R, (i, a, b))
        if not R.has(i): return R*eval_sum_symbolic(L, (i, a, b))
    if isinstance(f, C.Add):
        L, R = getab(f)
        lsum = eval_sum_symbolic(L, (i, a, b))
        rsum = eval_sum_symbolic(R, (i, a, b))
        if None not in (lsum, rsum):
            return lsum + rsum
    # Polynomial terms with Faulhaber's formula
    p = C.Wild('p')
    e = f.match(i**p)
    if e != None:
        c = p.subs_dict(e)
        B = C.bernoulli
        if c.is_integer and c >= 0:
            s = (B(c+1, b+1) - B(c+1, a))/(c+1)
            return s.expand()
    # Geometric terms
    if isinstance(f, C.Pow):
        r, k = f.args[:]
        if not r.has(i) and k == i:
            # TODO: Pow should be able to simplify x**oo depending
            # on whether |x| < 1 or |x| > 1 for non-rational x
            if (b is S.Infinity) and abs(r.evalf()) < 1:
                return r**a / (1-r)
            else:
                return (r**a - r**(b+1)) / (1-r)
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
