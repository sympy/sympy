"""1D quantum particle in a box."""

from sympy import Symbol, pi, sqrt, sin, Interval, S

from sympy.physics.quantum.operator import HermitianOperator
from sympy.physics.quantum.state import Ket, Bra, Wavefunction
from sympy.physics.quantum.constants import hbar
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.physics.quantum.hilbert import L2
from sympy.physics.quantum.cartesian import XKet

m = Symbol('m', real=True, positive=True)
L = Symbol('L', real=True, positive=True)


__all__ = [
    'PIABHamiltonian',
    'PIABKet',
    'PIABBra'
]


class PIABHamiltonian(HermitianOperator):
    """Particle in a box Hamiltonian operator."""

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    def _apply_operator_PIABKet(self, ket, **options):
        n = ket.label[0]
        return (n**2*pi**2*hbar**2)/(2*m*L**2)*ket


class PIABKet(Ket):
    """Particle in a box eigenket."""

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
        return PIABBra

    def _get_default_basis(self, **options):
        return XKet

    def _represent_XKet(self, basis, **options):
        from sympy.physics.quantum.represent import _append_index
        subs_info = options.get('subs',{})

        n = self.label[0]
        x = basis.position
        x = _append_index(x, options.pop("index", 1), **basis.label_assumptions)

        expr = sqrt(2/L)*sin(n*pi*x/L).subs(subs_info)
        return Wavefunction(expr, (x, 0, L))

    def _eval_innerproduct_PIABBra(self, bra):
        return KroneckerDelta(bra.label[0], self.label[0])


class PIABBra(Bra):
    """Particle in a box eigenbra."""

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    @classmethod
    def dual_class(self):
        return PIABKet

