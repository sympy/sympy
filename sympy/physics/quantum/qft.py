"""An implementation of qubits and gates acting on them."""

from sympy import Expr, Matrix, exp, I, pi
from sympy.matrices.matrices import eye
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import HilbertSpaceError
from sympy.physics.quantum.tensorproduct import matrix_tensor_product
from sympy.physics.quantum.applyops import apply_operators

from sympy.physics.quantum.gate import (
    Gate, HadamardGate, SwapGate, TwoQubitGate

)

__all__ = [
    'QFT',
    'IQFT',
]

#-----------------------------------------------------------------------------
# Fourier stuff
#-----------------------------------------------------------------------------


class RkGate(TwoQubitGate):
    """A Controlled phase gate.

    If Qubits specified in self.args[0] and self.args[1] are 1, then changes
    the phase of the state by e**(2*i*pi/2**k)

    *args are is the tuple describing which Qubits it should effect
    k is set by the third argument in the input, and describes how big of a
    phase shift it should apply

    >>> from sympy.physics.Qubit import Qubit, RkGate, apply_gates,\
    QubitZBasisSet
    >>> RkGate(1,0,2)
    R2(1, 0)
    >>> from sympy.physics.quantum import represent
    >>> represent(_, QubitZBasisSet(2))
    [1, 0, 0, 0]
    [0, 1, 0, 0]
    [0, 0, 1, 0]
    [0, 0, 0, I]
    >>> RkGate(1,0,3)*Qubit(1,1)
    R3(1, 0)*|11>
    >>> apply_gates(_)
    exp(pi*I/4)*|11>
    """
    gate_name = u'Rk'
    gate_name_latex = u'Rk'

    __slots__ = ['k']

    def __new__(cls, *args):
        obj = Gate.__new__(cls, *args[:-1])
        if 3 != len(args):
            num = obj.input_number
            raise QuantumError("This gate applies to %d Qubits" % (num))
        obj.k = args[-1]
        return obj

    def _apply_operator(self, Qubits):
        #switch Qubit basis and matrix basis when fully implemented
        mat = self.matrix
        args = [self.args[0][i] for i in reversed(range(2))]
        return self._apply(Qubits, mat, args)

    def _sympystr(self, printer, *args):
        return "R%s(%s, %s)" % (printer._print(self.k, *args),\
        printer._print(self.args[0][0], *args), printer._print(self.args[0][1], *args))

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,\
        exp(2*I*pi/2**self.k)]])

    @property
    def name(self):
        return "R%s(%s, %s)" % (self.k, self.args[0], self.args[1])

    @property
    def input_number(self):
        return 2

    def _print_operator_name_pretty(self, printer, *args):
        return prettyForm('R%s' % self.k)


class IRkGate(TwoQubitGate):
    """Inverse Controlled-Phase Gate

    Does the same thing as the RkGate, but rotates in the opposite direction
    within the complex plane. If Qubits specified in self.args[0]
    and self.args[1] are 1, then changes the phase of the state by
    e**(2*i*pi/2**k)

    *args are is the tuple describing which Qubits it should effect
    k is set by the third argument in the input, and describes how big of a
    phase shift it should apply

    >>> from sympy.physics.Qubit import Qubit, IRkGate, apply_gates,\
    QubitZBasisSet
    >>> IRkGate(1,0,2)
    IR2(1, 0)
    >>> from sympy.physics.quantum import represent
    >>> represent(_, QubitZBasisSet(2))
    [1, 0, 0,  0]
    [0, 1, 0,  0]
    [0, 0, 1,  0]
    [0, 0, 0, -I]
    >>> IRkGate(1,0,3)*Qubit(1,1)
    IR3(1, 0)*|11>
    >>> apply_gates(_)
    exp(-pi*I/4)*|11>
    """
    gate_name = u'IRk'
    gate_name_latex = u'IRk'

    def _sympystr(self, printer, *args):
        return "IR%s(%s, %s)" % (printer._print(self.k, *args),\
        printer._print(self.args[0][0], *args), printer._print(self.args[0][1], *args))

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,\
        exp(-2*I*pi/2**self.k)]])

    @property
    def name(self):
        return "IR%s(%s, %s)" % (self.k, self.args[0], self.args[1])

    def _print_operator_name_pretty(self, printer, *args):
        return prettyForm('IR%s' % self.k)


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

