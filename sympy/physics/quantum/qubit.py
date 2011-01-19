"""Qubits for quantum computing.

Todo:
* Implement measure in a new file, measure.py.
* Get IntQubit subclass working and integrate with Qubit.
* Update docstrings.
* Update tests.
"""

from sympy import Integer, log, Mul, Add, Pow
from sympy.core.basic import sympify
from sympy.core.containers import Tuple
from sympy.matrices.matrices import Matrix

from sympy.physics.quantum.hilbert import ComplexSpace
from sympy.physics.quantum.state import Ket, Bra, State

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.matrixcache import (
    numpy_ndarray, scipy_sparse_matrix
)

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
    def _eval_args(cls, args):
        # Turn strings into tuple of strings
        if len(args) == 1 and isinstance(args[0], basestring):
            args = tuple(args[0])

        args = sympify(args)

        # Validate input (must have 0 or 1 input)
        for element in args:
            if not (element == 1 or element == 0):
                raise ValueError("Qubit values must be 0 or 1, got: %r" % element)
        return args

    @classmethod
    def _eval_hilbert_space(cls, args):
        return ComplexSpace(2)**len(args)

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
        newargs = list(self.qubit_values)
        for i in bits:
            bit = int(self.dimension-i-1)
            if newargs[bit] == 1:
                newargs[bit] = 0
            else:
                newargs[bit] = 1
        return self.__class__(*tuple(newargs))


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
        """Represent this qubits in the computational basis (ZGate).
        """
        format = options.get('format', 'sympy')
        n = 1
        definite_state = 0
        for it in reversed(self.qubit_values):
            definite_state += n*it
            n = n*2
        result = [0]*(2**self.dimension)
        result[int(definite_state)] = 1
        if format == 'sympy':
            return Matrix(result)
        elif format == 'numpy':
            import numpy as np
            return np.matrix(result, dtype='complex').transpose()
        elif format == 'scipy.sparse':
            from scipy import sparse
            return sparse.csr_matrix(result, dtype='complex').transpose()


class QubitBra(QubitState, Bra):

    @property
    def dual_class(self):
        return Qubit

    def _represent_ZGate(self, basis, **options):
        format = options.get('format', 'sympy')
        result = self.dual._represent_ZGate(basis, **options)
        if format == 'sympy':
            return result.H
        elif format == 'numpy' or format == 'scipy.sparse':
            return result.transpose().conjugate()
        raise QuantumError('Invalide format: %r' % format)


def matrix_to_qubit(matrix):
    """Convert from the matrix repr. to a sum of Qubit objects.

    Parameters
    ----------
    matrix : Matrix
        The matrix to build the Qubit representation of.

    Examples
    --------
    """
    # Determine the format based on the type of the input matrix
    format = 'sympy'
    if isinstance(matrix, numpy_ndarray):
        format = 'numpy'
    if isinstance(matrix, scipy_sparse_matrix):
        format = 'scipy.sparse'

    # Make sure it is of correct dimensions for a Qubit-matrix representation.
    # This logic should work with sympy, numpy or scipy.sparse matrices.
    if matrix.shape[0] == 1:
        mlistlen = matrix.shape[1]
        nqubits = log(mlistlen, 2)
        ket = False
        cls = QubitBra
    elif matrix.shape[1] == 1:
        mlistlen = matrix.shape[0]
        nqubits = log(mlistlen, 2)
        ket = True
        cls = Qubit
    else:
        raise QuantumError(
            'Matrix must be a row/column vector, got %r' % matrix
        )
    if not isinstance(nqubits, Integer):
        raise QuantumError('Matrix must be a row/column vector of size '
                           '2**nqubits, got: %r' % matrix)
    # Go through each item in matrix, if element is non-zero, make it into a
    # Qubit item times the element.
    result = 0
    for i in range(mlistlen):
        if ket:
            element = matrix[i,0]
        else:
            element = matrix[0,i]
        if format == 'numpy' or format == 'scipy.sparse':
            element = complex(element)
        if element != 0.0:
            # Form Qubit array; 0 in bit-locations where i is 0, 1 in
            # bit-locations where i is 1
            qubit_array = [1 if i&(1<<x) else 0 for x in range(nqubits)]
            qubit_array.reverse()
            result = result + element*cls(*qubit_array)

    # If sympy simplified by pulling out a constant coefficient, undo that.
    if isinstance(result, (Mul,Add,Pow)):
        result = result.expand()

    return result


def qubit_to_matrix(qubit, format='sympy'):
    """Coverts an Add/Mul of Qubit objects into it's matrix representation

    This function is the inverse of ``matrix_to_qubit`` and is a shorthand
    for ``represent(qubit, ZGate(0))``.
    """
    from sympy.physics.quantum.gate import ZGate
    return represent(qubit, ZGate(0), format=format)

