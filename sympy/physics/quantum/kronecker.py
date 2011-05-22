"""Discrete Kronecker delta function."""

from sympy import S, Function
from sympy.printing.pretty.stringpict import prettyForm, stringPict

__all__ = [
    'KroneckerDelta'
]

class KroneckerDelta(Function):
    """The discrete, or Kronecker, delta function.

    A function that takes in two integers i and j. It returns 0 if i and j are
    not equal or it returns 1 if i and j are equal.

    Parameters
    ==========
    i : Number, Symbol
        The first index of the delta function.
    j : Number, Symbol
        The second index of the delta function.

    Examples
    ========

    A simple example with integer indices::

        >>> from sympy.physics.quantum import KroneckerDelta
        >>> KroneckerDelta(1,2)
        0
        >>> KroneckerDelta(3,3)
        1

    Symbolic indices::

        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> KroneckerDelta(i, j)
        d(i,j)
        >>> KroneckerDelta(i, i)
        1
        >>> KroneckerDelta(i, i+1)
        0
        >>> KroneckerDelta(i, i+1+k)
        d(i,i + k + 1)

    References
    ==========

    http://en.wikipedia.org/wiki/Kronecker_delta
    """

    nargs = 2
    is_commutative=True

    @classmethod
    def eval(cls, i, j):
        """
        Evaluates the discrete delta function.
        """
        if i > j:
            return cls(j,i)
        diff = i-j
        if diff == 0:
            return S.One
        elif diff.is_number:
            return S.Zero

    def _eval_subs(self, old, new):
        r = KroneckerDelta(self.args[0].subs(old, new), self.args[1].subs(old,\
        new))
        return r

    def _eval_dagger(self):
        return self

    def _latex_(self,printer):
        return "\\delta_{%s%s}"% (self.args[0].name,self.args[1].name)

    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)"% (self.__class__.__name__, self.args[0],\
        self.args[1])

    def _sympystr(self, printer, *args):
        return 'd(%s,%s)'% (self.args[0],self.args[1])

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right((prettyForm(','))))
        pform = prettyForm(*pform.right((printer._print(self.args[1], *args))))
        a = stringPict(u'\u03b4')
        b = pform
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _latex(self, printer, *args):
        i = printer._print(self.args[0], *args)
        j = printer._print(self.args[1], *args)
        return '\\delta_{%s %s}' % (i,j)
