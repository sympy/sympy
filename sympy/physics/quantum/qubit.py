"""Qubits for quantum computing.

Todo:
* Get things here to use a new MultiZGate as a basis for represent.
* Fix represent, matrix_to_qubit, qubit_to_matrix.
* Update docstrings.
"""

from sympy import Integer, I, log, Mul, Add, Pow
from sympy.core.basic import sympify
from sympy.core.containers import Tuple
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices.matrices import Matrix

from sympy.physics.quantum.hilbert import ComplexSpace
from sympy.physics.quantum.state import Ket, Bra, State

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.represent import represent

__all__ = [
    'Qubit',
    'Qubit',
    'QubitBra',
    'qubit_to_matrix',
    'matrix_to_qubit'
]

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
        """Returns the values of the qubits as a tuple."""
        return self.label

    #-------------------------------------------------------------------------
    # Special methods
    #-------------------------------------------------------------------------

    def __len__(self):
        return self.dimension

    def __getitem__(self, bit):
        return self.qubit_values[int(self.dimension-bit-1)]

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


class Qubit(QubitState, Ket):

    @property
    def dual_class(self):
        return QubitBra

    def _eval_innerproduct_QubitBra(self, bra, **hints):
        if self.label == bra.label:
            return Integer(1)
        else:
            return Integer(0)

    def _represent_ZGate(self, basis, **options):
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
        return Qubit


def matrix_to_qubit(matrix):
    """Converts a matrix representation of the state of a system into a Sum
    of Qubit objects

    Takes a matrix representation of a Qubit and puts in into a sum of Qubit
    eigenstates of the ZBasisSet. Can be used in conjunction with represent
    to turn the matrix representation returned into dirac notation.

    matrix argument is the matrix that shall be converted to the Add/Mul
    of Qubit eigenkets

    >>> from sympy.physics.Qubit import matrix_to_qubit, Qubit, QubitZBasisSet
    >>> from sympy.physics.quantum import represent
    >>> represent(Qubit(0,1), QubitZBasisSet(2))
    [0]
    [1]
    [0]
    [0]
    >>> matrix_to_qubit(represent(Qubit(0,1), QubitZBasisSet(2)))
    |'01'>
    """
    #make sure it is of correct dimensions for a Qubit-matrix representation
    qubit_number = log(matrix.rows,2)
    if matrix.cols != 1 or not isinstance(qubit_number, Integer):
        raise QuantumError()

    #go through each item in matrix, if element is not zero, make it into a
    #Qubit item times coefficient
    result = 0
    mlistlen = len(matrix.tolist())
    for i in range(mlistlen):
        if matrix[i] != 0:
            #form Qubit array; 0 in bit-locations where i is 0, 1 in
            #bit-locations where i is 1
            qubit_array = [1 if i&(1<<x) else 0 for x in range(qubit_number)]
            qubit_array.reverse()
            result = result + matrix[i]*Qubit(qubit_array)

    #if sympy simplified by pulling out a constant coefficient, undo that
    if isinstance(result, (Mul,Add,Pow)):
        result = result.expand()
    return result


def qubit_to_matrix(qubit):
    """Coverts a Add/Mul of Qubit objects into it's matrix representation

    This function takes in a dirac notation-like expression and express it
    as a Sympy matrix.

    >>> from sympy.physics.Qubit import qubit_to_matrix, Qubit
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> qubit_to_matrix(Qubit(0,0)/sqrt(2) + Qubit(0,1)/sqrt(2))
    [2**(1/2)/2]
    [2**(1/2)/2]
    [         0]
    [         0]
    """
    from sympy.physics.quantum.gate import ZGate
    return represent(qubit, ZGate(0))
