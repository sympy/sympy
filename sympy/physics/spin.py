from sympy import I, Symbol, S, Integer, Rational
from sympy.printing.pretty.stringpict import prettyForm, stringPict

from sympy.physics.quantum import HermitianOperator, State, Ket, Bra
from sympy.physics.quantum import KroneckerDelta, hbar
from sympy.physics.hilbert import ComplexSpace

class SpinBase(HermitianOperator):

    @classmethod
    def _eval_hilbert_space(cls, label):
        return ComplexSpace(2*label[0]+1)

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (self.__class__.__name__, printer._print(self.label,\
        *args))

    def _print_operator_name(self, printer, *args):
        return self.__class__.__name__

    def _print_operator_name_pretty(self, printer, *args):
        a = stringPict('S')
        b = stringPict(self._coord)
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _print_contents(self, printer, *args):
        return self._print_operator_name(printer, *args)

    def _print_contents_pretty(self, printer, *args):
        pform = self._print_operator_name_pretty(printer, *args)
        return pform


class Sx(SpinBase):

    _coord = 'x'

    def _eval_commutator(self, other):
        if isinstance(other, Sy):
            return I*hbar*Sz(self.label[0])
        elif isinstance(other, Sz):
            return -I*hbar*Sy(self.label[0])


class Sy(SpinBase):

    _coord = 'y'

    def _eval_commutator(self, other):
        if isinstance(other, Sx):
            return -I*hbar*Sz(self.label[0])
        elif isinstance(other, Sz):
            return I*hbar*Sy(self.label[0])


class Sz(SpinBase):

    _coord = 'z'

    def _eval_commutator(self, other):
        if isinstance(other, Sx):
            return -I*hbar*Sz(self.label[0])
        elif isinstance(other, Sz):
            return I*hbar*Sy(self.label[0])

    def _apply_to_ket_SzKet(self, ket):
        return (hbar*ket.label[1])*ket


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


class SpinState(State):
    pass


class SzKet(SpinState, Ket):

    _label_separator = ','

    def _apply_operator_Sz(self, op):
        return (hbar*self.label[1])*self

    def _apply_operator_S2(self, op):
        s = self.label[0]
        return (hbar*hbar*s*(s+1))

    def _eval_innerproduct(self, bra, **hints):
        if isinstance(bra, SzBra):
            d1 = KroneckerDelta(self.label[0], bra.label[0])
            d2 = KroneckerDelta(self.label[1], bra.label[1])
            return d1*d2

    @property
    def dual_class(self):
        return SzBra


class SzBra(SpinState, Bra):

    _label_separator = ','

    @property
    def dual_class(self):
        return SzKet
