from sympy import I, Symbol, S, Integer, Rational, Matrix, sqrt
from sympy.printing.pretty.stringpict import prettyForm, stringPict

from sympy.physics.quantum import HermitianOperator, State, Ket, Bra
from sympy.physics.quantum import KroneckerDelta, hbar
from sympy.physics.hilbert import ComplexSpace, HilbertSpaceError
from sympy.physics.matrices import msigma


class SpinBase(object):

    @property
    def s(self):
        return self.label[0]

    @classmethod
    def _eval_hilbert_space(cls, label):
        return ComplexSpace(2*label[0]+1)


class SpinOpBase(HermitianOperator, SpinBase):

    def _sympyrepr(self, printer, *args):
        return '%s(%s)' % (
            self.__class__.__name__, printer._print(self.label,*args)
        )

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

    def _check_hilbert_space_and_spin(self, other):
        if not self.hilbert_space == other.hilbert_space:
            raise HilbertSpaceError(
                'Incompatible hilbert spaces: %r and %r' % (self, other)
            )
        if not self.label[0] == other.label[0]:
            raise ValueError(
                'Spin value must be equal: %r != %r' % (
                    self.label[0], other.label[0]
                )
            )


class Sx(SpinOpBase):

    _coord = 'x'

    def _eval_commutator(self, other):
        if isinstance(other, Sy):
            return I*hbar*Sz(self.label[0])
        elif isinstance(other, Sz):
            return -I*hbar*Sy(self.label[0])

    def _represent_Sz(self, basis, **options):
        self._check_hilbert_space_and_spin(basis)

        if self.label[0] == Rational(1,2):
            return (hbar/2)*msigma(1)

    def _represent_Sx(self, basis, **options):
        self._check_hilbert_space_and_spin(basis)

        if self.label[0] == Rational(1,2):
            return (hbar/2)*msigma(3)


class Sy(SpinOpBase):

    _coord = 'y'

    def _eval_commutator(self, other):
        if isinstance(other, Sx):
            return -I*hbar*Sz(self.label[0])
        elif isinstance(other, Sz):
            return I*hbar*Sy(self.label[0])

    def _represent_Sz(self, basis, **options):
        self._check_hilbert_space_and_spin(basis)

        if self.label[0] == Rational(1,2):
            return (hbar/2)*msigma(2)

    def _represent_Sy(self, basis, **options):
        self._check_hilbert_space_and_spin(basis)

        if self.label[0] == Rational(1,2):
            return (hbar/2)*msigma(3)


class Sz(SpinOpBase):

    _coord = 'z'

    def _eval_commutator(self, other):
        if isinstance(other, Sx):
            return -I*hbar*Sz(self.label[0])
        elif isinstance(other, Sz):
            return I*hbar*Sy(self.label[0])

    def _apply_to_ket_SzKet(self, ket):
        return (hbar*ket.label[1])*ket

    def _represent_Sz(self, basis, **options):
        self._check_hilbert_space_and_spin(basis)

        if self.label[0] == Rational(1,2):
            return (hbar/2)*msigma(3)


class S2(SpinOpBase):

    def _eval_commutator(self, other):
        if isinstance(other, (Sx, Sy, Sz)):
            return S.Zero

    def _pretty(self, printer, *args):
        a = stringPict('S')
        b = stringPict('2')
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.above(top))


class SpinState(SpinBase, State):

    _label_separator = ','

    @property
    def ms(self):
        return self.label[1]


class SzKet(SpinState, Ket):

    def _apply_operator_Sz(self, op):
        return (hbar*self.ms)*self

    def _apply_operator_S2(self, op):
        s = self.s
        return (hbar*hbar*s*(s+1))

    def _eval_innerproduct(self, bra, **hints):
        if isinstance(bra, SzBra):
            d1 = KroneckerDelta(self.s, bra.s)
            d2 = KroneckerDelta(self.ms, bra.ms)
            return d1*d2

    @property
    def dual_class(self):
        return SzBra

    def _represent_Sz(self, basis, **options):
        if self.s == Rational(1,2):
            if self.ms == Rational(1,2):
                return Matrix([1,0])
            elif self.ms == -Rational(1,2):
                return Matrix([0,1])


class SzBra(SpinState, Bra):

    @property
    def dual_class(self):
        return SzKet

    def _represent_Sz(self, basis, **options):
        return self.dual._represent_Sz(basis, **options).H

class SxKet(SpinState, Ket):

    @property
    def dual_class(self):
        return SxBra

    def _apply_operator_S2(self, op):
        s = self.s
        return (hbar*hbar*s*(s+1))

    def _eval_innerproduct(self, bra, **hints):
        if isinstance(bra, SxBra):
            d1 = KroneckerDelta(self.s, bra.s)
            d2 = KroneckerDelta(self.ms, bra.ms)
            return d1*d2

    def _represent_Sz(self, basis, **options):
        if self.s == Rational(1,2):
            if self.ms == Rational(1,2):
                return Matrix([1,1])/sqrt(2)
            elif self.ms == -Rational(1,2):
                return Matrix([1,-1])/sqrt(2)


class SxBra(SpinState, Bra):

    @property
    def dual_class(self):
        return SxKet

    def _represent_Sz(self, basis, **options):
        return self.dual._represent_Sz(basis, **options).H

class SyKet(SpinState, Ket):

    @property
    def dual_class(self):
        return SyBra

    def _apply_operator_S2(self, op):
        s = self.s
        return (hbar*hbar*s*(s+1))

    def _eval_innerproduct(self, bra, **hints):
        if isinstance(bra, SyBra):
            d1 = KroneckerDelta(self.s, bra.s)
            d2 = KroneckerDelta(self.ms, bra.ms)
            return d1*d2

    def _represent_Sz(self, basis, **options):
        if self.s == Rational(1,2):
            if self.ms == Rational(1,2):
                return Matrix([1,I])/sqrt(2)
            elif self.ms == -Rational(1,2):
                return Matrix([1,-I])/sqrt(2)


class SyBra(SpinState, Bra):

    @property
    def dual_class(self):
        return SyKet

    def _represent_Sz(self, basis, **options):
        return self.dual._represent_Sz(basis, **options).H
