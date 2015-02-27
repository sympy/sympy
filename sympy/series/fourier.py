"""Fouries series expansion"""

from sympy import oo, pi
from sympy.core.compatibility import range
from sympy.core.sympify import sympify
from sympy.core.symbol import Symbol
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.compatibility import integer_types
from sympy.integrals import integrate
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.concrete.summations import Sum

class SeqFormula(object):
    """
    Generates coefficients of a sequence based on a formula

    Examples
    ========

    >>> from sympy.abc import x, n
    >>> from sympy import oo
    >>> from sympy.series.fourier import SeqFormula
    >>> s = SeqFormula(x*(x+1), x, (0, oo))
    >>> s.coeff(n)
    n*(n + 1)
    >>> s.coeff(10)
    110

    """
    def __init__(self, formula, x, bounds=None):
        self._formula = sympify(formula)
        self.x = sympify(x)
        self._coeff = {}
        if not bounds:
            bounds = (0, oo)
        self._lower, self._upper = bounds[0], bounds[1]
        self._current = bounds[0]

    def coeff(self, index):
        """
        returns the nth coefficient
        """
        if isinstance(index, integer_types):
            if index > self._upper or index < self._lower:
                raise IndexError("index %d not in (%d, %d)" % (index, self._lower, self._upper))
        if index in self._coeff:
            return self._coeff[index]
        self._coeff[index] = self._eval_coeff(index)
        return self._coeff[index]

    def __iter__(self):
        return self

    def _eval_coeff(self, index):
        return self._formula.subs(self.x, index)

    def next(self):
        if self._current > self._upper:
            raise StopIteration()
        val = self._eval_coeff(self._current)
        self._coeff[self._current] = val
        self._current += 1
        return val
    __next__ = next

class FourierSeries(Expr):
    """Represents Fourier sine/cosine series

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy import pi
    >>> from sympy.series.fourier import FourierSeries
    >>> fseries = FourierSeries(x, x, (-pi, pi))
    >>> fseries.as_series()
    Sum(-2*(-1)**n*sin(n*x)/n, (n, 1, oo))
    >>> fseries.as_nseries(n=3)
    2*sin(x) - sin(2*x) + 2*sin(3*x)/3

    """
    def __new__(cls, f, x, bounds=(-pi, pi)):
        f = sympify(f)
        x = sympify(x)

        lower_bound, upper_bound = bounds
        L = upper_bound - lower_bound
        int_tuple = (x, lower_bound, upper_bound)

        n = Symbol('n', integer=True, positive=True)
        a0 = 2 * (integrate(f, int_tuple) / L)

        fsubbed = f.subs(x, -x) # to check odd or even

        if fsubbed == -f:
            seq_an = SeqFormula(0, n, (1, oo)) # if odd
        else:
            cos_term = cos(2*n*pi*x / L)
            seq_an = SeqFormula((2*integrate(f * cos_term, int_tuple) / L), n, (1, oo))

        if fsubbed == f:
            seq_bn = SeqFormula(0, n, (1, oo)) # if even
        else:
            sin_term = sin(2*n*pi*x / L)
            seq_bn = SeqFormula((2*integrate(f * sin_term, int_tuple) / L), n, (1, oo))

        obj = Expr.__new__(cls, f, x)
        obj.x, obj.n, obj.bounds, obj.L = x, n, bounds, L
        obj.seq_an, obj.seq_bn = seq_an, seq_bn
        obj.a0, obj.an, obj.bn = a0, seq_an.coeff(n), seq_bn.coeff(n)
        return obj

    def free_symbols(self):
        return self.x.free_symbols

    def __getitem__(self, index):
        if isinstance(index, integer_types):
            if index < 0:
                raise IndexError("Index should be in (0, Inf)")
            return self._eval_n_term(index)
        elif isinstance(index, slice):
            return [self._eval_n_term(i) for i in range(index.start, index.stop, index.step or 1)]

    def as_series(self):
        """
        Represents series as:
        (a0 / 2) + Sum(an*cos(n*x), (n,1,oo)) + Sum(bn*sin(n*x), (n,1,oo))
        """
        n = self.n
        cos_term = cos(2*n*pi*self.x / self.L)
        sin_term = sin(2*n*pi*self.x / self.L)
        an = self.seq_an.coeff(n).simplify()
        bn = self.seq_bn.coeff(n).simplify()
        fseries = self.a0 / 2
        if an != 0:
            fseries += Sum(an * cos_term, (n, 1, oo))
        if bn != 0:
            fseries += Sum(bn * sin_term, (n, 1, oo))
        return fseries

    def _eval_n_term(self, n):
        """
        evaluates the nth term
        """
        cos_term = cos(2*n*pi*self.x / self.L)
        sin_term = sin(2*n*pi*self.x / self.L)
        if n == 0:
            return self.a0 / 2
        an = self.seq_an.coeff(n)
        bn = self.seq_bn.coeff(n)
        return Add(*[an*cos_term, bn*sin_term])

    def as_nseries(self, n=6):
        """
        Gives the terms till n = num
        """
        return Add(*self[0:n + 1])

def fourier_series(f, x, bounds=(-pi, pi), n=6):
    """
    Gives fourier series expansion
    if n is given, n terms will be returned
    if n is None, returns an infinite generator of terms
    """
    fseries = FourierSeries(f, x, bounds)
    if not n is None:
        return fseries.as_nseries(n)
    else:
        def yield_fourier_series(fseries):
            i = 0
            while True:
                yield fseries[i]
                i += 1
        return yield_fourier_series(fseries)
