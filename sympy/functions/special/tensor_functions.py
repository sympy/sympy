import operator
from sympy.core.compatibility import reduce
from sympy.core.function import Function
from sympy.core import sympify, S, Integer
from sympy.core.mul import prod

###############################################################################
###################### Kronecker Delta, Levi-Civita etc. ######################
###############################################################################

class Dij(Function):
    """
    Represents the Kronecker Delta Function

    if i == j, Dij(i, j) = 1
    otherwise Dij(i, j) = 0
    where i, j are usually integers
    """
    nargs = (1, 2)

    @classmethod
    def eval(cls, i, j=0):
        i, j = map(sympify, (i, j))
        if i == j:
            return S.One
        elif i.is_number and j.is_number:
            return S.Zero

def Eijk(*args, **kwargs):
    """
    Represent the Levi-Civita symbol.

    This is just compatibility wrapper to LeviCivita().
    """
    return LeviCivita(*args, **kwargs)

def eval_levicivita(*args):
    """Evaluate Levi-Civita symbol."""
    from sympy import factorial
    n = len(args)
    return prod(
        prod(args[j] - args[i] for j in xrange(i + 1, n))
        / factorial(i) for i in xrange(n))
    # converting factorial(i) to int is slightly faster

class LeviCivita(Function):
    """Represent the Levi-Civita symbol.

    For even permutations of indices it returns 1, for odd permutations -1, and
    for everything else (a repeated index) it returns 0.

    Thus it represents an alternating pseudotensor.

    >>> from sympy import LeviCivita, symbols
    >>> LeviCivita(1,2,3)
    1
    >>> LeviCivita(1,3,2)
    -1
    >>> LeviCivita(1,2,2)
    0
    >>> i,j,k = symbols('i j k')
    >>> LeviCivita(i,j,k)
    LeviCivita(i, j, k)
    >>> LeviCivita(i,j,i)
    0
    """
    @classmethod
    def eval(cls, *args):
        if all(isinstance(a, (int, Integer)) for a in args):
            return eval_levicivita(*args)
        if len(set(args)) < len(args):
            return S.Zero

    def doit(self):
        return eval_levicivita(*self.args)

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

