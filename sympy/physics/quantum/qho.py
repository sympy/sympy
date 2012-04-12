"""1D quantum harmonic oscillator."""

from sympy import Symbol, pi, sqrt, sin, Interval, S, factorial, exp
from sympy.functions import hermite

from sympy.physics.quantum.operator import HermitianOperator, Operator
from sympy.physics.quantum.state import Ket, Bra, Wavefunction
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.kronecker import KroneckerDelta
from sympy.physics.quantum.hilbert import L2
from sympy.physics.quantum.cartesian import XKet

m = Symbol('m', real=True, positive=True)
w = Symbol('w', real=True, positive=True) #omega


__all__ = [
    'QHOHamiltonian',
    'QHOKet',
    'QHOBra',
    'AOp',
    'ADaggerOp'
]


class QHOHamiltonian(HermitianOperator):
    """Quantum harmonic oscillator Hamiltonian operator."""

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    def _apply_operator_QHOKet(self, ket, **options):
        n = ket.label[0]
        return hbar*w*(n+S(1)/2)*ket


class QHOKet(Ket):
    """Quantum harmonic oscillator eigenket.

    Can be represented in continuous bases and manipulated with continuous
    operators as well.

    """

    @classmethod
    def default_args(self):
        return ('n',)

    @classmethod
    def def_label_assumptions(self):
        return {'integer':True}

    @classmethod
    def _eval_hilbert_space(cls, args):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    @classmethod
    def dual_class(self):
        return QHOBra

    def _get_default_basis(self, **options):
        return XKet

    def _represent_XKet(self, basis, **options):
        from sympy.physics.quantum.represent import _append_index
        subs_info = options.get('subs',{})

        n = self.label[0]
        x = basis.position
        x = _append_index(x, options.pop("index", 1), **basis.label_assumptions)

        prefactor = sqrt(S(1)/(2**n*factorial(n)))*(m*w/(pi*hbar))**(S(1)/4)
        expo = exp(-m*w*x**2/(2*hbar))
        herm = hermite(n, sqrt(m*w/hbar)*x)
        return Wavefunction(prefactor*expo*herm, x)

    def _eval_innerproduct_QHOBra(self, bra):
        return KroneckerDelta(bra.label[0], self.label[0])


class QHOBra(Bra):
    """Particle in a box eigenbra."""

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    @classmethod
    def dual_class(self):
        return QHOKet

class AOp(Operator):
    """ QHO Lowering operator """

    @classmethod
    def default_args(self):
        ('a',)

    def _apply_operator_QHOKet(self, ket, **options):
        n = ket.label[0]
        if n == 0:
            return 0

        return sqrt(n)*QHOKet(n-1)

    def _eval_dagger(self):
        return ADaggerOp(*self.args)

    def _eval_commutator_ADaggerOp(self, ket, **options):
        return S.One

class ADaggerOp(Operator):
    """ QHO raising operator """

    @classmethod
    def default_args(self):
        ('a_d',)

    def _apply_operator_QHOKet(self, ket, **options):
        n = ket.label[0]
        return sqrt(n+1)*QHOKet(n+1)

    def _eval_dagger(self):
        return AOp(*self.args)
