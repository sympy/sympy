from sympy import I, Symbol, S
from sympy.printing.pretty.stringpict import prettyForm, stringPict

from sympy.physics.quantum import Operator
from sympy.physics.hilbert import l2

hbar = Symbol('hbar')

class SpinBase(Operator):

    @classmethod
    def _eval_hilbert_space(cls, label):
        return l2(2*label+1)

    def _sympystr(self, printer, *args):
        return printer._print(self.__class__.__name__, *args)

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, printer._print(self.label,\
        *args))

    def _pretty(self, printer, *args):
        a = stringPict('S')
        b = stringPict(self._coord)
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))


class Sx(SpinBase):

    _coord = 'x'

    def _eval_commutator(self, other):
        if isinstance(other, Sy):
            return I*hbar*Sz(self.label)
        elif isinstance(other, Sz):
            return -I*hbar*Sy(self.label)


class Sy(SpinBase):

    _coord = 'y'

    def _eval_commutator(self, other):
        if isinstance(other, Sx):
            return -I*hbar*Sz(self.label)
        elif isinstance(other, Sz):
            return I*hbar*Sy(self.label)


class Sz(SpinBase):

    _coord = 'z'

    def _eval_commutator(self, other):
        if isinstance(other, Sx):
            return -I*hbar*Sz(self.label)
        elif isinstance(other, Sz):
            return I*hbar*Sy(self.label)

class S2(SpinBase):

    def _eval_commutator(self, other):
        if isinstance(other, (Sx, Sy, Sz)):
            return S.Zero

    def _pretty(self, printer, *args):
        a = stringPict('S')
        b = stringPict('2')
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.above(top))
