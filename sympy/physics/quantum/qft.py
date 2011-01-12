"""An implementation of qubits and gates acting on them."""

from sympy import Expr
from sympy.matrices.matrices import eye

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import HilbertSpaceError
from sympy.physics.quantum.tensorproduct import matrix_tensor_product
from sympy.physics.quantum.applyops import apply_operators

from sympy.physics.quantum.gate import (
    Gate, HadamardGate, RkGate, SwapGate, IRkGate

)

__all__ = [
    'QFT',
    'IQFT',
]

#-----------------------------------------------------------------------------
# Fourier stuff
#-----------------------------------------------------------------------------

class Fourier(Gate):
    """Superclass of Quantum Fourier and Inverse Quantum Fourier Gates

    This gate represents the quantum fourier tranform. It can be decomposed
    into elementary gates using the famous QFT decomposition.

    Takes in two args telling which Qubits to start and stop doing the (I)QFT

    >>> from sympy.physics.Qubit import Qubit, QFT, IQFT, QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> represent(QFT(0,2).decompose(), QubitZBasisSet(2))
    [1/2,  1/2,  1/2,  1/2]
    [1/2,  I/2, -1/2, -I/2]
    [1/2, -1/2,  1/2, -1/2]
    [1/2, -I/2, -1/2,  I/2]
    >>> QFT(0,2).decompose()
    SwapGate(0,1)*HadamardGate(0)*R2(1, 0)*HadamardGate(1)
    >>> represent(IQFT(0,2).decompose(), QubitZBasisSet(2))
    [1/2,  1/2,  1/2,  1/2]
    [1/2, -I/2, -1/2,  I/2]
    [1/2, -1/2,  1/2, -1/2]
    [1/2,  I/2, -1/2, -I/2]
    >>> IQFT(0,2).decompose()
    HadamardGate(1)*IR2(1, 0)*HadamardGate(0)*SwapGate(0,1)
    """

    def __new__(self, *args):
        if args[0] >= args[1]:
            raise QuantumError("Start must be smaller than finish")
        return Expr.__new__(self, *args)

    @property
    def _apply(self, Qubits):
        raise NotImplementedError("This command shouldn't happen")

    @property
    def minimum_dimension(self):
        #Can apply to a Qubit basisSet up to one less that its last arg
        return self.args[1]-1

    @property
    def input_number(self):
        #first input should be start of register, second is end of register
        return 2

    def _represent_ZBasisSet(self, hilbert_size, format = 'sympy'):
        if hilbert_size <= self.minimum_dimension:
            raise HilbertSpaceError("hilbert_size doesn't work")
        product = []
        product.append(eye(2**(self.args[0])))
        product.append(self.matrix)
        product.append(eye(2**(hilbert_size - self.args[1])))
        return matrix_tensor_product(*product)

    @property
    def matrix(self):
        # Can't form the matrix on its own yet
        NotImplementedError("Fourier Transforms don't know how yet")

    def _apply_ZBasisSet(self, Qubits):
        #decomposes self into constituients and applies
        return apply_operators(self.decompose*Qubits)


class QFT(Fourier):

    def decompose(self):
        """Decomposes QFT into elementary gates."""
        start = self.args[0]
        finish = self.args[1]
        circuit = 1
        for level in reversed(range(start, finish)):
            circuit = HadamardGate(level)*circuit
            for i in range(level-start):
                circuit = RkGate(level, level-i-1, i+2)*circuit
        for i in range((finish-start)/2):
            circuit = SwapGate(i+start, finish-i-1)*circuit
        return circuit


class IQFT(Fourier):

    def decompose(self):
        """Decomposes IQFT into elementary gates."""
        start = self.args[0]
        finish = self.args[1]
        circuit = 1
        for i in range((finish-start)/2):
            circuit = SwapGate(i+start, finish-i-1)*circuit
        for level in range(start, finish):
            for i in reversed(range(level-start)):
                circuit = IRkGate(level, level-i-1, i+2)*circuit
            circuit = HadamardGate(level)*circuit
        return circuit

