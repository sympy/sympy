"""An implementation of qubits and gates acting on them.

Todo:
* Update docstrings.
* Update tests.
* Implement apply using decompose.
* Implement represent using decompose or something smarter. For this to work
  we first have to implement represent for SWAP.
"""

from sympy import Expr, Matrix, exp, I, pi, Integer
from sympy.matrices.matrices import eye
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import HilbertSpaceError
from sympy.physics.quantum.tensorproduct import matrix_tensor_product
from sympy.physics.quantum.applyops import apply_operators

from sympy.physics.quantum.gate import (
    Gate, HadamardGate, SwapGate, OneQubitGate, CGate, PhaseGate, TGate, ZGate
)

__all__ = [
    'QFT',
    'IQFT',
    'RkGate',
    'Rk'
]

#-----------------------------------------------------------------------------
# Fourier stuff
#-----------------------------------------------------------------------------


class RkGate(OneQubitGate):
    """This is the R_k gate of the QTF."""

    gate_name = u'Rk'
    gate_name_latex = u'R'

    def __new__(cls, *args, **old_assumptions):
        if len(args) != 2:
            raise QuantumError(
                'Rk gates only take two arguments, got: %r' % args
            )
        # For small k, Rk gates simplify to other gates, using these 
        # substitutions give us familiar results for the QFT for small numbers
        # of qubits.
        target = args[0]
        k = args[1]
        if k == 1:
            return ZGate(target)
        elif k == 2:
            return PhaseGate(target)
        elif k == 3:
            return TGate(target)
        args = cls._eval_args(args)
        inst = Expr.__new__(cls, *args, **{'commutative':False})
        # Now set the slots on the instance
        inst.hilbert_space = cls._eval_hilbert_space(args)
        return inst

    @property
    def k(self):
        return self.label[1]

    @property
    def targets(self):
        return self.label[:1]

    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[1,0],[0,exp(Integer(2)*pi*I/(Integer(2)**self.k))]])
        raise NotImplementedError('Invalid format for the R_k gate: %r' % format)


Rk = RkGate


class Fourier(Gate):
    """Superclass of Quantum Fourier and Inverse Quantum Fourier Gates."""

    @classmethod
    def _eval_args(self, args):
        if len(args) != 2:
            raise QuantumError(
                'QFT/IQFT only takes two arguments, got: %r' % args
            )
        if args[0] >= args[1]:
            raise QuantumError("Start must be smaller than finish")
        return Gate._eval_args(args)

    @property
    def targets(self):
        return range(self.label[0],self.label[1])


class QFT(Fourier):
    """The forward quantum Fourier transform."""

    gate_name = u'QFT'
    gate_name_latex = u'QFT'

    def decompose(self):
        """Decomposes QFT into elementary gates."""
        start = self.label[0]
        finish = self.label[1]
        circuit = 1
        for level in reversed(range(start, finish)):
            circuit = HadamardGate(level)*circuit
            for i in range(level-start):
                circuit = CGate(level-i-1, RkGate(level, i+2))*circuit
        for i in range((finish-start)/2):
            circuit = SwapGate(i+start, finish-i-1)*circuit
        return circuit

    def _eval_inverse(self):
        return IQFT(*self.args)


class IQFT(Fourier):
    """The inverse quantum Fourier transform."""

    gate_name = u'IQFT'
    gate_name_latex = u'{QFT^{-1}}'

    def decompose(self):
        """Decomposes IQFT into elementary gates."""
        start = self.args[0]
        finish = self.args[1]
        circuit = 1
        for i in range((finish-start)/2):
            circuit = SwapGate(i+start, finish-i-1)*circuit
        for level in range(start, finish):
            for i in reversed(range(level-start)):
                circuit = CGate(level-i-1, RkGate(level, -i-2))*circuit
            circuit = HadamardGate(level)*circuit
        return circuit

    def _eval_inverse(self):
        return QFT(*self.args)
