__all__ = ['DynamicSymbol']

from sympy import Symbol, S

class DynamicSymbol(Symbol):
    """A symbol implicitly dependent upon time.

    DynamicSymbol behaves just like Symbol, except when it is differentitated
    with respect to Symbol('t'), in which case  another DynamicSymbol is
    returned. The new DynamicSymbol has a 'd' appended to the end of its name,
    to represent "time derivative", with one 'd' per level of differentiation.
    To avoid confusion, it is reccommended that DynamicSymbols are not given
    base names which end in 'd'.

    >>> t, x = symbols('t x')
    >>> y = DynamicSymbol('y')
    >>> diff(2*x + 3*y, t)
    3*yd
    >>> isinstance(diff(y, t), DynamicSymbol)
    True

    Additionally, DynamicSymbol enables you to take the derivative of a Sympy
    expression with respect to a DynamicSymbol.

    >>> diff(2*x**2 + 3*y, x)
    4*x

    This is in contrast to differentiating a Sympy Function, which
    returns a Derivative object. A Sympy expression can't be differentiated
    with respect to a Derivative object.

    This functionality is essential for implementing Lagrange's method,
    Hamilton's method, and Kane's method, all of which rely on differentiating
    with respect to time varying quantities.  For example, in Lagrange's
    method, you must form :math:`\frac{\partial L}{\partial\dot{q}}`, which
    won't work in Sympy if :math:`\dot{q}}` is represented as a Derivative
    object.

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

if __name__ == "__main__":
    import doctest
    from sympy import symbols, diff
    global_dict = {'symbols': symbols,
                   'diff' : diff,
                   'DynamicSymbol': DynamicSymbol}
    doctest.testmod(globs=global_dict)
