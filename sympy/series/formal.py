from __future__ import print_function, division

from sympy import C, S, symbols, Symbol, collect, Function, Add
from sympy.core.sympify import sympify
from sympy.core.relational import Eq
from sympy.functions.elementary.piecewise import Piecewise
from sympy.concrete import Sum
from sympy.core import Add, Pow
from sympy.core.expr import Expr
from sympy.polys.partfrac import apart
from sympy.solvers import solve


class FormalSeries(object):
    def __new__(self, *args, **kwargs):
        formula = kwargs.pop("formula", None)
        if formula:
            return FormalSeriesFormula(args[0], formula)

        sequence = kwargs.pop('sequence', None)
        if sequence:
            return FormalSeriesSeq(args[0], sequence)

        function = kwargs.pop('function', None)
        if function:
            return FormalSeriesFunc(args[0], function)


class FormalSeriesBase(Expr):
    """
    Formal Power Series base class
    """
    def __add__(self, other):
        it = self.it
        gen = self.gen + (other.gen).subs(other.it, it)
        return FormalSeriesFormula(self.sym, (gen, it))

    def __sub__(self, other):
        it = self.it
        gen = self.gen - (other.gen).subs(other.it, it)
        return FormalSeriesFormula(self.sym, (gen, it))

    def __getitem__(self, key):
        gen, it = self.gen, self.it
        sym = self.sym
        if isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or 6
            step = key.step or 1
            return [self[k] for k in range(start, stop, step)]
        else:
            return gen.subs(it, key) * sym**key

    def as_series(self, n=6):
        sym = self.sym
        s = self[0:n]
        return Add(*s) + C.Order(sym**n)


class FormalSeriesFunc(FormalSeriesBase):
    def __init__(self, sym, function):
        self.gen, self.it = self.findgen(sym, function)
        self.sym = sym

    def findgen(self, sym, function):
        for k in range(0, 5):
            fdiff = function.diff(sym, k)
            if fdiff.is_rational_function():
                return self.findgenR(sym, function)  # Integrate k times

        DE, f = self.simpleDE(sym, function)
        RE = self.DEtoRE(sym, DE, f)
        k = Symbol('k', integer=True)
        return Piecewise((S.One/k, True)), k

    def simpleDE(self, sym, function):
        a = symbols('a:4')
        f = Function('f')(sym)

        # Search for case when k=1
        eq, DE = self.makeDE(sym, function, f, a, S.One)
        sol = solve(eq, a[:1], dict=True)
        if sol and sol[0].values()[0].is_rational_function():
            return DE.subs(sol[0].keys()[0], sol[0].values()[0])

        for k in range(2, 5):
            eq, DE = self.makeDE(sym, function, f, a, k)
            terms = [t.as_independent(sym)[1] for t in eq.as_ordered_terms()]
            eq = collect(eq, [function.diff(sym, k) for k in range(0, k)])
            coeff = [t.as_independent(sym)[0] for t in eq.as_ordered_terms()]
            if self.independentK(sym, terms, k):
                sol = solve(coeff, a, dict=True)[0]
                for key, value in sol.iteritems():
                    DE = DE.subs(key, value)
                return DE, f

        raise NotImplementedError('Cannot find simple DE')

    def independentK(self, sym, terms, k):
        summands = []
        for s in terms:
            summands += s.as_ordered_terms()
        summands = list(set(summands))
        l = len(summands)
        for i in range(l):
            for j in range(i+1, l):
                if (summands[i] / summands[j]).is_rational_function():
                    return False
        return True

    def makeDE(self, sym, function, f, a, order):
        eq = function.diff(sym, order)
        DE = f.diff(sym, order)
        for k in range(0, order):
            eq += a[k] * function.diff(sym, k)
            DE += a[k] * f.diff(sym, k)
        return eq, DE

    def DEtoRE(self, sym, DE, f):
        print (DE)
        return None

    def findgenR(self, sym, function):
        gen = S.Zero
        k = Symbol('k', integer=True)
        terms = apart(function).as_ordered_terms()
        for t in terms:
            c, d = t.as_numer_denom()
            d, j = d.as_base_exp()
            a = -d.as_coeff_add()[0]
            gen += (-1)**j * c * C.binomial(j+k-1, k).rewrite(C.factorial) / a**(j+k)
        return Piecewise((gen, True)), k


class FormalSeriesFormula(FormalSeriesBase):
    def __init__(self, sym, formula):
        self.gen, self.it = self.findgen(sym, formula)
        self.sym = sym

    def findgen(self, sym, formula):
        if type(formula) != tuple:
            raise TypeError("Formula should be 'tuple'")
        gen, it = formula
        if gen.is_Piecewise:
            return gen, it
        return Piecewise((gen, True)), it


class FormalSeriesSeq(FormalSeriesBase):
    def __init__(self, sym, sequence):
        self.gen, self.it = self.findgen(sym, sequence)
        self.sym = sym

    def findgen(self, sym, sequence):
        k = Symbol('k', integer=True)
        l = len(sequence)
        cond = [(sequence[i], Eq(k%l, i)) for i in range(l)]
        return Piecewise(*cond), k
