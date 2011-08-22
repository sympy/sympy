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
    """Particle in a box eigenket.

    Can be represented in continuous bases and manipulated with continuous
    operators as well.

    Examples
    ========

    >>> from sympy.physics.quantum.piab import PIABKet, PIABBra
    >>> from sympy.physics.quantum.cartesian import XOp, XKet, PxOp

    A default instance will simply give the ket an index of n, where n is an
    integer

    >>> p = PIABKet()
    >>> p
    |n>
    >>> p.label[0].is_integer
    True

    When you represent the ket, you get a Wavefunction object. By default, it
    will be represented in the XKet (1D cartesian) basis

    >>> from sympy.physics.quantum import represent
    >>> represent(p)
    Wavefunction(2**(1/2)*sin(pi*n*x/L)/L**(1/2), (x, 0, L))

    You can also then compute more complicated representations.
    In this example, we get an extra factor of x in the resulting Wavefunction,
    as expected.

    >>> represent(XOp()*p, basis=XKet)
    Wavefunction(2**(1/2)*x*sin(pi*n*x/L)/L**(1/2), (x, 0, L))

    In this example, the momentum operator in the position basis is a
    differential operator, so we actually end up with the derivative of the
    original represented expression.

    >>> represent(PxOp()*p, basis=XKet)
    Wavefunction(-2**(1/2)*hbar*I*pi*n*cos(pi*n*x/L)/L**(3/2), (x, 0, L))

    We can even compute expectation values for this state!

    >>> p_bra = PIABBra()
    >>> represent(p_bra*XOp()*p, basis=XKet)
    L*cos(pi*n)**2/2

    cos(pi*n)**2 is 1, so this gives L/2, as expected!

    >>> represent(p_bra*PxOp()*p, basis=XKet)
    0

    The expectation value of the momentum for a piab state is 0, as expected!
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

