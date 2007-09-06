
from sympy.core import Basic, S, Add, Mul, Symbol, Equality, Interval
from sympy.core.methods import NoRelMeths, ArithMeths

class Sum(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated summation."""

    precedence = Basic.Apply_precedence

    def __new__(cls, f, *symbols, **assumptions):
        f = Basic.sympify(f)

        if isinstance(f, Basic.Number):
            if isinstance(f, Basic.NaN):
                return S.NaN
            elif isinstance(f, Basic.Zero):
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
                            limits.append(tuple(V))
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
        if not hints.get('sums', True):
            return self

    def _eval_summation(self, f, x):
        return

def sum(*args, **kwargs):
    summation = Sum(*args, **kwargs)

    if isinstance(summation, Sum):
        return summation.doit()
    else:
        return summation
