"""An implementation of qubits and gates acting on them."""

from sympy import Expr, Add, Mul, Pow, Integer, I, pi
from sympy.core.numbers import Number
from sympy.core.basic import sympify
from sympy.core.containers import Tuple
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices.matrices import Matrix, eye
from sympy.printing.pretty.stringpict import prettyForm, stringPict

from sympy.physics.qexpr import QuantumError, QExpr
from sympy.physics.hilbert import ComplexSpace, HilbertSpaceError
from sympy.physics.quantum import (
    Operator, Ket, Bra, State,
    represent, apply_operators,
    matrix_tensor_product
)

#-----------------------------------------------------------------------------
# Qubit Classes
#-----------------------------------------------------------------------------




class QubitState(State):
    """Represents a single definite eigenstate of the ZBasisSet

    Qubit object contains a tuple of values which represent the quantum
    numbers of each Qubit. It can have individual bits flipped such that a
    one becomes a zero and vice versa. This Object is also callable allowing
    user to pick out the value of a particular bit. The Qubit class can also
    have outDecimal class variable set to True which causes the state to
    present itself in decimal form.

    Object can be instantiated by different args:
        - *args with each representing the value of a single Qubit
        - A single decimal value. Code will convert to binary using least
          number of bits possible
        - A decimal value and the number of bits you wish to express it in
          (The size of the Hilbert Space)

    >>> from sympy.physics.Qubit import Qubit
    >>> Qubit(0,0,0)
    |'000'>
    >>> Qubit(5)
    |'101'>
    >>> a = Qubit(5,4)
    >>> a
    |'0101'>
    >>> a.flip(0)
    |'0100'>
    >>> len(a)
    4
    >>> a.dimension
    4
    >>> a[0]
    1
    >>> a.name
    '5'
    """

    #-------------------------------------------------------------------------
    # Initialization/creation
    #-------------------------------------------------------------------------
    
    @classmethod
    def _eval_label(cls, label):
        # Turn string into tuple of ints, sympify. This must be a tuple to
        # make qubits hashable.
        label = sympify(tuple(label))
        
        # Validate input (must have 0 or 1 input)
        if isinstance(label, (str, list, tuple)):
            for element in label:
                if not (element == 1 or element == 0):
                    raise ValueError("Qubit values must be 0 or 1, got: %r" % element)
        else:
            raise TypeError('str, list, tuple expected, got: %r' % label)
        return Tuple(*label)

    @classmethod
    def _eval_hilbert_space(cls, label):
        return ComplexSpace(2)**len(label)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def dimension(self):
        """The number of Qubits in the state."""
        return len(self.qubit_values)

    @property
    def nqubits(self):
        return self.dimension

    @property
    def qubit_values(self):
        """Returns the values of the qubits as a list."""
        return self.label

    @property
    def x_basis_transform(self):
        #Transform matrix from ZBasis to XBasis and back again
        return 1/sqrt(2)*Matrix([[1,1],[1,-1]])

    @property
    def y_basis_transform(self):
        #Transform matrix from ZBasis to YBasis and back again
        return Matrix([[I,0],[0,-I]])

    #-------------------------------------------------------------------------
    # Special methods
    #-------------------------------------------------------------------------

    def __len__(self):
        return self.dimension

    def __getitem__(self, bit):
        if bit > self.dimension - 1:
            raise IndexError('bit index out of range: %r' % bit)
        return self.self.qubit_values[int(self.dimension-bit-1)]

    #-------------------------------------------------------------------------
    # Printing methods
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # Utility methods
    #-------------------------------------------------------------------------

    def flip(self, *bits):
        """Flip the bit(s) given."""
        newargs = list(self.label)
        for i in bits:
            bit = int(self.dimension-i-1)
            if newargs[bit] == 1:
                newargs[bit] = 0
            else:
                newargs[bit] = 1
        return self.__class__(newargs)


class QubitKet(QubitState, Ket):

    @property
    def dual_class(self):
        return QubitBra

    def _eval_innerproduct_QubitBra(self, bra, **hints):
        if self.label == bra.label:
            return Integer(1)
        else:
            return Integer(0)

    def _represent_MultiZGate(self, basis, **options):
        """
        if basis.hilbert_space != l2(2)**self.dimension:
            raise HilbertSpaceException("Basis and Qubit\
            dimensions do not match!")
        TODO Figure out validation here
        """
        n = 1
        definiteState = 0
        for it in reversed(self.qubit_values):
            definiteState += n*it
            n = n*2
        result = [0]*(2**self.dimension)
        result[int(definiteState)] = 1
        return Matrix(result)


class QubitBra(QubitState, Bra):

    @property
    def dual_class(self):
        return QubitKet


# Shorthand notation for the ket.
Qubit = QubitKet

#-----------------------------------------------------------------------------
# Gate Super-Classes
#-----------------------------------------------------------------------------

class Gate(Operator):
    """A unitary operator that acts on qubits.

    The superclass of gate operators that acts on qubit(s)
        - minimum_dimension is the least size a set fo Qubits must be for
          the gate to be applied
        - input number is the number of Qubits the gate acts on (e.g. a
          three Qubit gate acts on three Qubits)
        - name returns the string of the Qubit
        - apply applies gate onto Qubits
        - _represent_Qubit?BasisSet represents the gates matrix in the basis
    """

    _label_separator = ','

    gate_name = u'G'
    gate_name_pretty = u'G'
    gate_name_latex = u'G'

    nqubits = Integer(1)
    __slots__ = ['min_qubits']

    #-------------------------------------------------------------------------
    # Initialization/creation
    #-------------------------------------------------------------------------

    def __new__(cls, label, **assumptions):
        gate = Operator.__new__(cls, label, **assumptions)
        gate.min_qubits = max(gate.label)
        return gate

    @classmethod
    def _eval_label(cls, label):
        label = Operator._eval_label(label)
        temp_nqubits = len(label)
        if not temp_nqubits == cls.nqubits:
            raise QuantumError(
                'Gate acts on %r qubits, but got %r qubit indices' % \
                (cls.nqubits, temp_nqubits)
            )
        for i in range(temp_nqubits):
            if not label[i].is_Integer:
                raise TypeError('Integer expected, got: %r' % label[i])
            if (label[i] in label[:i]) or (label[i] in label[i+1:]):
                raise QuantumError('Qubit indices in a gate cannot be duplicated')
        return label

    @classmethod
    def _eval_hilbert_space(cls, label):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(label)+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def matrix(self):
        raise NotImplementedError("matrixRep Not implemented")

    @property
    def x_basis_transform(self):
        #Transform matrix from ZBasis to XBasis and back again
        return 1/sqrt(2)*Matrix([[1,1],[1,-1]])

    @property
    def y_basis_transform(self):
        #Transform matrix from ZBasis to YBasis and back again
        return Matrix([[I,0],[0,-I]])

    #-------------------------------------------------------------------------
    # Print methods
    #-------------------------------------------------------------------------

    def _print_contents(self, printer, *args):
        label = self._print_label(printer, *args)
        return '%s(%s)' % (self.gate_name, label)

    def _print_contents_pretty(self, printer, *args):
        a = stringPict(unicode(self.gate_name_pretty))
        b = self._print_label_pretty(printer, *args)
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _latex(self, printer, *args):
        label = self._print_label(printer, *args)
        return '%s_{%s}' % (self.gate_name_latex, label)

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    def _apply_operator_QubitKet(self, qubits, **options):
        """Apply this gate to a QubitKet."""
        mat = self.matrix
        label = self.label

        # Check number of qubits this gate acts on.
        if qubits.nqubits < self.min_qubits:
            raise HilbertSpaceError(
                'Gate needs a minimum of %r qubits to act on, got: %r' %\
                    (self.min_qubits, qubits.nqubits)
            )

        # Find which column of the matrix this Qubit applies to.
        column_index = 0
        n = 1
        for element in label:
            column_index += n*qubits[element]
            n = n<<1
        column = mat[:,int(column_index)]

        # Now apply each column element to Qubit.
        result = 0
        for index in range(len(column.tolist())):
            new_qubit = Qubit(*qubits.args)
            # Flip the bits that need to be flipped.
            for bit in range(len(label)):
                if new_qubit[label[bit]] != (index>>bit)&1:
                    new_qubit = new_qubit.flip(label[bit])
            # The value in that row and column times the flipped-bit qubit
            # is the result for that part.
            result += column[index]*new_qubit
        return result

    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    def _represent(self, basis, format='sympy'):
        if isinstance(basis, Gate):
            basis_size = 1
        elif isinstance(basis, Pow) and isinstance(basis.base, Gate):
            basis_size = basis.exp
        else:
            raise QuantumError(
                'Basis must be a gate operator, or tensor product of gate operators.'
            )

        if basis_size < self.min_qubits:
            raise QuantumError(
                'The basis given is too small for the given Gate objects.'
            )

        gate = self.matrix
        if isinstance(basis, Gate):
            return gate
        else:
            m = represent_hilbert_space(
                gate, basis_size, self.label, format
            )
            return m


#-----------------------------------------------------------------------------
# Single Qubit Gates
#-----------------------------------------------------------------------------

class OneQubitGate(Gate):

    nqubits = Integer(1)


class HadamardGate(OneQubitGate):
    """An object representing a Hadamard Gate

    This gate puts creates an even superposition of eigenstates.
    It maps the state |0> -> |0>/sqrt(2) + |1>/sqrt(2)
    and |1> -> |0>/sqrt(2) - |1>/sqrt(2). Can be applied or represented
    using apply_gates and represent.

    >>> from sympy.physics.Qubit import HadamardGate, Qubit, apply_gates,
    >>> HadamardGate(0)*Qubit(0,0)
    HadamardGate(0)*|00>
    >>> apply_gates(_)
    2**(1/2)/2*|00> + 2**(1/2)/2*|01>
    >>> from sympy.physics.quantum import represent
    >>> represent(HadamardGate(0)*Qubit(0,0), QubitZBasisSet(2))
    [2**(1/2)/2]
    [2**(1/2)/2]
    [         0]
    [         0]
    """
    gate_name = u'H'
    gate_name_pretty = u'H'
    gate_name_latex = u'H'
        
    @property
    def matrix(self):
        from sympy.functions.elementary.miscellaneous import sqrt
        return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))


class XGate(OneQubitGate):
    """An object representing a Pauli-X gate (AKA NOT Gate)

    This is a NOT Gate because it flips the value of a bit to its opposite
    just like the classical-NOT gate. Thus it maps the states:
    |1> -> |0> and |0> -> |1>


    >>> from sympy.physics.Qubit import Qubit, XGate, apply_gates,\
    QubitZBasisSet
    >>> XGate(0)*Qubit(0,0)
    XGate(0)*|00>
    >>> apply_gates(_)
    |01>
    >>> from sympy.physics.quantum import represent
    >>>represent(XGate(0)*Qubit(0,0), QubitZBasisSet(2))
    >>> represent(XGate(0)*Qubit(0,0), QubitZBasisSet(2))
    [1]
    [1]
    [0]
    [0]
    """
    gate_name = u'X'
    gate_name_pretty = u'X'
    gate_name_latex = u'X'
            
    @property
    def matrix(self):
        return Matrix([[0, 1], [1, 0]])


class YGate(OneQubitGate):
    """An object representing a Pauli-Y gate

    This gate flips the bit given in self.args[0] and then applies a relative
    phase shift (pi/2 if 1, and -pi/2 if 0).
    Thus it maps the states:
    |0> -> I*|1> and |1> -> -I*|0>

    >>> from sympy.physics.Qubit import Qubit, YGate, apply_gates, \
    QubitZBasisSet
    >>> YGate(0)*Qubit(0,1)
    YGate(0)*|01>
    >>> apply_gates(_)
    -1.0*I*|00>
    >>> from sympy.physics.quantum import represent
    >>> represent(YGate(0)*Qubit(0,0), QubitZBasisSet(2))
    [0]
    [I]
    [0]
    [0]
    """
    gate_name = u'Y'
    gate_name_pretty = u'Y'
    gate_name_latex = u'Y'
        
    @property
    def matrix(self):
        return Matrix([[0, complex(0,-1)], [complex(0,1), 0]])


class ZGate(OneQubitGate):
    """An object representing a Pauli-Z gate

    This gate rotates the relative phase of an eigenstate by pi if the state
    is 1. Thus the gate maps the states:
    |0> -> |0> and |1> -> -|1>

    >>> from sympy.physics.Qubit import Qubit, ZGate, apply_gates, \
    QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> ZGate(0)*Qubit(0,1)
    ZGate(0)*|01>
    >>> apply_gates(_)
    -1*|01>
    >>> represent(ZGate(0)*Qubit(0,0), QubitZBasisSet(2))
    [1]
    [0]
    [0]
    [0]
    """
    gate_name = u'Z'
    gate_name_pretty = u'Z'
    gate_name_latex = u'Z'

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, -1]])


class PhaseGate(OneQubitGate):
    """An object representing a phase gate

    This gate rotates the phase of the eigenstate by pi/2 if the state is 1.
    Thus the gate maps the states:
    |0> -> |0> and |1> -> I*|1>

    >>> from sympy.physics.Qubit import Qubit, PhaseGate, apply_gates,\
    QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> represent(PhaseGate(0)*Qubit(0,1), QubitZBasisSet(2))
    [0]
    [I]
    [0]
    [0]
    >>> apply_gates(PhaseGate(0)*Qubit(0,1))
    I*|01>
    >>> PhaseGate(0)*Qubit(0,1)
    PhaseGate(0)*|01>
    """
    gate_name = u'S'
    gate_name_pretty = u'S'
    gate_name_latex = u'S'

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, complex(0,1)]])


class TGate(OneQubitGate):
    """An object representing a pi/8 gate

    This gate rotates the phase of the eigenstate by pi/2 if the state is 1.
    Thus the gate maps the states:
    |0> -> |0> and |1> -> exp(I*pi/4)*|1>

    >>> from sympy.physics.Qubit import Qubit, TGate, apply_gates,\
    QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> TGate(0)*Qubit(0,1)
    TGate(0)*|01>
    >>> apply_gates(_)
    exp(pi*I/4)*|01>
    >>> represent(TGate(0)*Qubit(0,1), QubitZBasisSet(2))
    [          0]
    [exp(pi*I/4)]
    [          0]
    [          0]
    """
    gate_name = u'T'
    gate_name_pretty = u'T'
    gate_name_latex = u'T'

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, exp(I*pi/4)]])

#-----------------------------------------------------------------------------
# 2 Qubit Gates
#-----------------------------------------------------------------------------

class TwoQubitGate(Gate):

    nqubits = Integer(2)


class RkGate(Gate):
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
    gate_name_pretty = u'Rk'
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


class IRkGate(RkGate):
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
    gate_name_pretty = u'IRk'
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


class CTGate(Gate):
    """Controlled Pi/8 Gate (Controlled Version of TGate)

    Applies a TGate if the Qubit specified by self.args[0] is True.
    Thus, it rotates the phase by pi/2 if both Qubits are true.
    It maps the state:
    |11> -> exp(I*pi/4)*|11> (leaves others unaffected)

    >>> from sympy.physics.Qubit import Qubit, CTGate, apply_gates,\
    QubitZBasisSet
    >>> apply_gates(CTGate(0,1)*Qubit(1,1))
    exp(pi*I/4)*|11>
    """
    gate_name = u'CT'
    gate_name_pretty = u'CT'
    gate_name_latex = u'CT'

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,\
        exp(I*pi/4)]])


class CZGate(Gate):
    """Controlled Z-Gate

    Applies a ZGate if the Qubit specified by self.args[0] is true. Thus, it
    rotates the phase by pi if both Qubits are true.
    i.e It maps the state: |11> -> -|11>

    >>> from sympy.physics.Qubit import Qubit, CZGate, apply_gates,\
    QubitZBasisSet
    >>> apply_gates(CZGate(0,1)*Qubit(1,1))
    -1*|11>
    """
    gate_name = u'CZ'
    gate_name_pretty = u'CZ'
    gate_name_latex = u'CZ'
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])


class CPhaseGate(Gate):
    """Controlled Phase-Gate

    Applies the phase gate contingent on the value specified in self.args[0]
    being true. Thus, it rotates the phase by pi/2 if both Qubits are true.
    i.e. It maps the state: |11> -> exp(I*pi/4)*|11>

    >>> from sympy.physics.Qubit import Qubit, CPhaseGate, apply_gates,\
    QubitZBasisSet
    >>> apply_gates(CPhaseGate(0,1)*Qubit(1,1))
    I*|11>
    """
    gate_name = 'CS'
    gate_name_pretty = u'CS'
    gate_name_latex = u'CS'

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,I]])


class CNotGate(Gate):
    """Controlled NOT-Gate

    Note: (This is the 'entangling gate' most often made reference to in the
    literature) CNot, Hadamard and the Pauli-Gates make
    a universal group in quantum computation. Flips the second Qubit
    (The target) contingent on the first Qubit being 1. Can be thought of as
    a reversible XOR Gate.

    >>> from sympy.physics.Qubit import Qubit, CNotGate, apply_gates,\
    QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> represent(CNotGate(0,1), QubitZBasisSet(2))
    [1, 0, 0, 0]
    [0, 0, 0, 1]
    [0, 0, 1, 0]
    [0, 1, 0, 0]
    >>> apply_gates(CNotGate(0,1)*Qubit(0,1))
    |'11'>
    """
    gate_name = 'CNot'
    gate_name_pretty = u'CNot'
    gate_name_latex = u'CNot'

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])


class SwapGate(Gate):
    """A SwapGate Object

    Thus swaps two Qubits locations within the Tensor Product.

    >>> from sympy.physics.Qubit import Qubit, SwapGate, apply_gates,\
    QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> represent(SwapGate(0,1), QubitZBasisSet(2))
    [1, 0, 0, 0]
    [0, 0, 1, 0]
    [0, 1, 0, 0]
    [0, 0, 0, 1]
    >>> apply_gates(SwapGate(0,1)*Qubit(0,1))
    |'10'>
    """
    gate_name = 'Swap'
    gate_name_pretty = u'Swap'
    gate_name_latex = u'Swap'

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])


class ToffoliGate(Gate):
    """A ToffoliGate (AKA the Controlled-Controlled (double controlled) NOT

    Flips the third Qubit (the target) contingent on the first and
    second Qubits being 1 (first and second Qubits are control bits)
    It can be thought of as a Controlled CNotGate

    >>> from sympy.physics.Qubit import Qubit, ToffoliGate, apply_gates,\
    QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> represent(ToffoliGate(1,0,2), QubitZBasisSet(3))
    [1, 0, 0, 0, 0, 0, 0, 0]
    [0, 1, 0, 0, 0, 0, 0, 0]
    [0, 0, 1, 0, 0, 0, 0, 0]
    [0, 0, 0, 0, 0, 0, 0, 1]
    [0, 0, 0, 0, 1, 0, 0, 0]
    [0, 0, 0, 0, 0, 1, 0, 0]
    [0, 0, 0, 0, 0, 0, 1, 0]
    [0, 0, 0, 1, 0, 0, 0, 0]
    >>> apply_gates(ToffoliGate(1,0,2)*Qubit(1,1,1))
    |'011'>
    """
    gate_name = u'Toffoli'
    gate_name_pretty = u'Toffoli'
    gate_name_latex = u'Toffoli'

    @property
    def matrix(self):
        return Matrix([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],\
        [0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1]\
        ,[0,0,0,0,0,0,1,0]])

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

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------


def _operator(first, second, format):
    """Returns the Outer product of a one or zero ket and bra"""
    if (first != 1 and first != 0) or (second != 1 and second != 0):
        raise QuantumError("can only make matricies |0><0|, |1><1|, |0><1|,\
        or |1><0|")
    if first:
        if second:
            ret = Matrix([[0,0],[0,1]])
        else:
            ret = Matrix([[0,0],[1,0]])
    else:
        if second:
            ret = Matrix([[0,1],[0,0]])
        else:
            ret = Matrix([[1,0],[0,0]])
    if format == 'sympy':
        return ret
    else:
        import numpy as np
        return np.matrix(ret.tolist())


def _np_tensor_product(*product):
    """
        Wrapper Function that abstracts away numpy kron function as a
        tensor_product function
    """
    import numpy as np
    answer = product[0]
    for item in product[1:]:
        answer = np.kron(answer, item)
    return answer


def represent_hilbert_space(gateMatrix, hilbert_size, Qubits, format='sympy'):
    """
        This is a helper function used by Gates represent functions to represent
        their matricies in a given HilbertSpace (i.e. a Hadamard Gate matrix
        looks different if applied to a different number of Qubits)

        Inputs: gateMatrix      -- The fundamental input matrix
                hilbert_size    -- the number of qbits in the hilbertspace
                Qubits          -- the qubit(s) we will be applying to 
    """
    if format == 'sympy':
        eye = getattr(gateMatrix, 'eye')
        kron = matrix_tensor_product
    elif format=='numpy':
        #if user specified numpy as matrix format, try to import
        try:
            import numpy as np
        except Exception:
            #If we couldn't load, just revert to sympy
            represent_hilbert_space(gateMatrix, hilbert_size,\
            Qubits, format='sympy')
        #redirect eye to np.eye function, and kron to a modified numpy function
        gateMatrix = np.matrix(gateMatrix.tolist())
        aeye = getattr(np, 'eye')
        eye = lambda x: np.matrix(aeye(x))
        kron = _np_tensor_product
    else:
        raise ValueError()

    if gateMatrix.shape[1] == 2:
        product = []
        Qubit = Qubits[0]
        #fill product with [I1,Gate,I2] such that the unitaries,
        # I, cause the gate to be applied to the correct Qubit
        if Qubit != hilbert_size-1:
            product.append(eye(2**(hilbert_size-Qubit-1)))
        product.append(gateMatrix)
        if Qubit != 0:
            product.append(eye(2**Qubit))

        #do the tensor product of these I's and gates
        if format == 'sympy' or format == 'numpy':
            MatrixRep = kron(*product)
        else:
            raise ValueError()
        return MatrixRep

    #If we are dealing with a matrix that is inheritely multi-qubit
    else:
        #find the control and target Qubit(s)
        controls = Qubits[:-1]
        controls = [x for x in reversed(controls)]
        target =  Qubits[-1]
        answer = 0
        product = []
        #break up gateMatrix into list of 2x2 matricies's
        #This list will be used for determining what matrix goes where
        matrixArray = []
        for i in range(gateMatrix.shape[1]/2):
            for j in range(gateMatrix.shape[1]/2):
                matrixArray.append(gateMatrix[i*2:i*2+2,j*2:j*2+2])

        #Build up tensor products and additions, so that we can form matrix
        for i in range((gateMatrix.shape[1]/2)**2):
            product = []
            #Put Unities in all locations
            for j in range(hilbert_size):
                product.append(eye(2))
            n = 0
            #put Operators |0><0|, |1><1|, |0><1|, or |1><0|
            # in place of I's for all control bits
            for item in controls:
                product.pop(hilbert_size-1-item)
                #Operator is picked so that if i = 0xyab
                #(base 2; x,y,a,b = 0 or 1),
                #then operator is |xy><ab| = |x><a|X|y><b| each of (|x><a|)
                # which goes in control bit location
                product.insert(hilbert_size-1-item,\
                 _operator(i>>(n+len(controls))&1,(i>>n)&1, format))
                n = n+1
            #put the correct submatrix from matrixarray into target-bit location
            product.pop(hilbert_size-1-target)
            product.insert(hilbert_size-1-target, matrixArray[i])

            #preform Tensor product first time
            if isinstance(answer, (int, Integer)):
                if format == 'sympy' or format == 'numpy':
                    answer = kron(*product)
                else:
                    raise ValueError()
            #add last answer to tensor_product of what we have
            else:
                if format == 'sympy' or format == 'numpy':
                    answer = answer + kron(*product)
                else:
                    raise ValueError()
        return answer


def matrix_to_Qubits(matrix):
    """
        Converts a matrix representation of the state of a system into a Sum
        of Qubit objects

        Takes a matrix representation of a Qubit and puts in into a sum of Qubit
        eigenstates of the ZBasisSet. Can be used in conjunction with represent
        to turn the matrix representation returned into dirac notation.

        matrix argument is the matrix that shall be converted to the Add/Mul
        of Qubit eigenkets

        >>> from sympy.physics.Qubit import matrix_to_Qubits, Qubit, QubitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(Qubit(0,1), QubitZBasisSet(2))
        [0]
        [1]
        [0]
        [0]
        >>> matrix_to_Qubits(represent(Qubit(0,1), QubitZBasisSet(2)))
        |'01'>
    """
    #make sure it is of correct dimensions for a Qubit-matrix representation
    Qubit_number = log(matrix.rows,2)
    if matrix.cols != 1 or not isinstance(Qubit_number, Integer):
        raise QuantumError()

    #go through each item in matrix, if element is not zero, make it into a
    #Qubit item times coefficient
    result = 0
    mlistlen = len(matrix.tolist())
    for i in range(mlistlen):
        if matrix[i] != 0:
            #form Qubit array; 0 in bit-locations where i is 0, 1 in
            #bit-locations where i is 1
            Qubit_array = [1 if i&(1<<x) else 0 for x in range(Qubit_number)]
            Qubit_array.reverse()
            result = result + matrix[i]*Qubit(Qubit_array)

    #if sympy simplified by pulling out a constant coefficient, undo that
    if isinstance(result, (Mul,Add,Pow)):
        result = result.expand()
    return result


def Qubits_to_matrix(Qubits):
    """
        Coverts a Add/Mul of Qubit objects into it's matrix representation

        This function takes in a dirac notation-like expression and express it
        as a Sympy matrix.

        >>> from sympy.physics.Qubit import Qubits_to_matrix, Qubit
        >>> from sympy.functions.elementary.miscellaneous import sqrt
        >>> Qubits_to_matrix(Qubit(0,0)/sqrt(2) + Qubit(0,1)/sqrt(2))
        [2**(1/2)/2]
        [2**(1/2)/2]
        [         0]
        [         0]
    """
    return represent(Qubits, ZGate(0))


def gate_simp(circuit):
    """Simplifies gates symbolically

    It first sorts gates using gate_sort. It then applies basic
    simplification rules to the circuit, e.g., XGate**2 = Identity

    >>> from sympy.physics.Qubit import gate_simp, HadamardGate
    >>> gate_simp(HadamardGate(1)**3*HadamardGate(0))
    HadamardGate(0)*HadamardGate(1)
    """

    #bubble sort out gates that commute
    circuit = gate_sort(circuit)

    #do simplifications by subing a simplification into the first element
    #which can be simplified
    #We recursively call gate_simp with new circuit as input
    #more simplifications exist
    if isinstance(circuit, Mul):
        #Iterate through each element in circuit; simplify if possible
        for i in range(len(circuit.args)):
            #H,X,Y or Z squared is 1. T**2 = S, S**2 = Z
            if isinstance(circuit.args[i], Pow):
                if isinstance(circuit.args[i].base, \
                    (HadamardGate, XGate, YGate, ZGate))\
                    and isinstance(circuit.args[i].exp, Number):
                    #Build a new circuit taking replacing the 
                    #H, X,Y,Z squared with one
                    newargs = (circuit.args[:i] + (circuit.args[i].base**\
                    (circuit.args[i].exp % 2),) + circuit.args[i+1:])
                    #Recursively simplify the new circuit
                    circuit = gate_simp(Mul(*newargs))
                    break
                elif isinstance(circuit.args[i].base, PhaseGate):
                    #Build a new circuit taking old circuit but splicing 
                    #in simplification
                    newargs = circuit.args[:i]
                    #replace PhaseGate**2 with ZGate
                    newargs = newargs + (ZGate(circuit.args[i].base.args[0][0])**\
                    (Integer(circuit.args[i].exp/2)), circuit.args[i].base**\
                    (circuit.args[i].exp % 2))
                    #append the last elements
                    newargs = newargs + circuit.args[i+1:]
                    #Recursively simplify the new circuit
                    circuit =  gate_simp(Mul(*newargs))
                    break
                elif isinstance(circuit.args[i].base, TGate):
                    #Build a new circuit taking all the old elements
                    newargs = circuit.args[:i]
                    
                    #put an Phasegate in place of any TGate**2
                    newargs = newargs + (SGate(circuit.args[i].base.args[0][0])**\
                    Integer(circuit.args[i].exp/2), circuit.args[i].base**\
                    (circuit.args[i].exp % 2))
                    
                    #append the last elements
                    newargs = newargs + circuit.args[i+1:]
                    #Recursively simplify the new circuit
                    circuit =  gate_simp(Mul(*newargs))
                    break

    return circuit


def gate_sort(circuit):
    """Sorts the gates while keeping track of commutation relations

    This function uses a bubble sort to rearrange the order of gate
    application. Keeps track of Quantum computations special commutation
    relations (e.g. things that apply to the same Qubit do not commute with
    each other)

    circuit is the Mul of gates that are to be sorted.

    >>> from sympy.physics.Qubit import HadamardGate, XGate, YGate, CNotGate,\
    gate_sort
    >>> gate_sort(YGate(2)**2*HadamardGate(0)*CNotGate(0,1)*XGate(1)*YGate(0))
    HadamardGate(0)*CNotGate(0,1)*YGate(0)*XGate(1)*(YGate(2))**2
    """
    #bubble sort of gates checking for commutivity of neighbor
    changes = True
    while changes:
        changes = False
        cirArray = circuit.args
        for i in range(len(cirArray)-1):
            #Go through each element and switch ones that are in wrong order
            if isinstance(cirArray[i], (Gate, Pow)) and\
            isinstance(cirArray[i+1], (Gate, Pow)):
                #If we have a Pow object, look at only the base
                if isinstance(cirArray[i], Pow):
                    first = cirArray[i].base
                else:
                    first = cirArray[i]

                if isinstance(cirArray[i+1], Pow):
                    second = cirArray[i+1].base
                else:
                    second = cirArray[i+1]

                #If the elements should sort
                if first.args[0][0] > second.args[0][0]:
                    #make sure elements commute
                    #meaning they do not affect ANY of the same Qubits
                    commute = True
                    for arg1 in first.args[0]:
                       for arg2 in second.args[0]:
                            if arg1 == arg2:
                                commute = False
                    # if they do commute, switch them
                    if commute:
                        circuit = Mul(*(circuit.args[:i] + (circuit.args[i+1],)\
                         + (circuit.args[i],) + circuit.args[i+2:]))
                        cirArray = circuit.args
                        changes = True
                        break
    return circuit
