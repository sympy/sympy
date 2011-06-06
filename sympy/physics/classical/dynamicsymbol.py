__all__ = ['DynamicSymbol']

from sympy import Symbol, S

class DynamicSymbol(Symbol):
    """
    Class for time-varying quantities.  When DynamicSymbol's derivative
    is taken with respect to Symbol 't', a time differentiated version is
    returned.
    """

    @property
    def free_symbols(self):
        return set([Symbol('t'), self])

    def _eval_derivative(self, s):
        if s == Symbol('t'):
            return DynamicSymbol(self.name + 'd')
        elif self == s:
            return S.One
        else:
            return S.Zero
