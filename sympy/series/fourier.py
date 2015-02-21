"""Fouries series expansion"""

from sympy.core import oo, Symbol

class SeqFormula(object):
    """
    Creates a generator of coefficients based on a formula
    """
    def __init__(self, formula, x, bounds=None):
        self._formula = formula
        self.x = x
        self._coeff = {}
        if not bounds:
            bounds = (0, oo)
        self._lower, self._upper = bounds[0], bounds[1]
        self._current = bounds[0]

    def coeff(self, index):
        if index > self._upper or index < self._lower:
            raise IndexError("index %d not in (%d, %d)" % (index, self._lower, self._upper))
        if self._coeff.has_key(index):
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
