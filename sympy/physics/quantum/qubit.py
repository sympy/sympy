"""Qubits for quantum computing.

Todo:
* Finish implementing measurement logic. This should include partial
  measurements as well as POVM.
* Update docstrings.
* Update tests.
"""

import math

from sympy import Integer, log, Mul, Add, Pow, conjugate
from sympy.core.basic import sympify
from sympy.core.containers import Tuple
from sympy.matrices.matrices import Matrix, zeros
from sympy.printing.pretty.stringpict import prettyForm
from sympy.functions.elementary.miscellaneous import sqrt


from sympy.physics.quantum.hilbert import ComplexSpace
from sympy.physics.quantum.state import Ket, Bra, State

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.matrixcache import (
    numpy_ndarray, scipy_sparse_matrix
)

__all__ = [
    'Qubit',
    'QubitBra',
    'IntQubit',
    'IntQubitBra',
    'qubit_to_matrix',
    'matrix_to_qubit',
    'measure_all',
    'measure_partial',
    'measure_partial_oneshot',
    'measure_all_oneshot'
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
        # If we are passed a QubitState or subclass, we just take its qubit
        # values directly.
        if len(args) == 1 and isinstance(args[0], QubitState):
            return args[0].qubit_values

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


class IntQubitState(QubitState):
    """A base class for qubits that work with binary representations of ints."""

    @classmethod
    def _eval_args(cls, args):
        # The case of a QubitState instance
        if len(args) == 1 and isinstance(args[0], QubitState):
            return QubitState._eval_args(args)
        # For a single argument, we construct the binary representation of
        # that integer with the minimal number of bits.
        if len(args) == 1 and args[0] > 1:
            #rvalues is the minimum number of bits needed to express the number
            rvalues = reversed(
                range(int(math.ceil(math.log(args[0], 2)+.01)+.001))
            )
            qubit_values = [(args[0]>>i)&1 for i in rvalues]
            return QubitState._eval_args(qubit_values)
        # For two numbers, the second number is the number of bits
        # on which it is expressed, so IntQubit(0,5) == |00000>.
        elif len(args) == 2 and args[1] > 1:
            #TODO Raise error if there are not enough bits
            qubit_values = [(args[0]>>i)&1 for i in reversed(range(args[1]))]
            return QubitState._eval_args(qubit_values)
        else:
            return QubitState._eval_args(args)

    def as_int(self):
        number = 0
        n = 1
        for i in reversed(self.qubit_values):
            number += n*i
            n = n<<1
        return number

    def _print_label(self, printer, *args):
        return str(self.as_int())

    def _print_label_pretty(self, printer, *args):
        label = self._print_label(printer, *args)
        return prettyForm(label)

    _print_label_repr = _print_label
    _print_label_latex = _print_label


class IntQubit(IntQubitState, Qubit):

    @property
    def dual_class(self):
        return IntQubitBra


class IntQubitBra(IntQubitState, QubitBra):

    @property
    def dual_class(self):
        return IntQubit


#-----------------------------------------------------------------------------
# Qubit <---> Matrix conversion functions
#-----------------------------------------------------------------------------


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


#-----------------------------------------------------------------------------
# Measurement
#-----------------------------------------------------------------------------


def measure_all(qubit, format='sympy'):
    """Given a qubit, return the primitive states and their probabilities."""
    m = qubit_to_matrix(qubit, format)
    
    if format == 'sympy':
        results = []
        m = m.normalized()
        size = max(m.shape) #max of shape to account for bra or ket
        nqubits = int(math.log(size)/math.log(2))
        for i in range(size):
            if m[i] != 0.0:
                results.append(
                    (Qubit(IntQubit(i, nqubits)), m[i]*conjugate(m[i]))
                )
        return results
    else:
        raise NotImplementedError("This function can't handle non-sympy" +\
                                  "matrix formats yet")   
                                  
def measure_partial(qubit, bit, format='sympy'):
    """Does a partial ensemble measure on the specifed qubit 
       TODO Make work for multiple input bits
    """
    m = qubit_to_matrix(qubit, format)

    if format == 'sympy':
        result = []
        m = m.normalized()
        size = max(m.shape) #max of shape to account for bra or ket
        nqubits = int(math.log(size,2)+.1)
        bit_mask = 1<<bit #when bit_mask is '&' with an array index,
        #It will determine if the bit is true in that state
        
        #break the matrix into halves one where specified bit is true, other not
        true_matrix  = zeros((2**nqubits, 1))
        false_matrix = zeros((2**nqubits, 1)) 
        for i in range(2**nqubits):
            if bit_mask&i:
                true_matrix[i] = m[i]
            else:
                false_matrix[i] = m[i]

        #calculate probability of finding the specified bit as true
        prob_true = 0
        for item in true_matrix*true_matrix.H:
            prob_true += item
        prob_false = 1 - prob_true
        
        #If the qubit observed is already in a definite state, return input
        if prob_true == 0 or prob_false == 0:
            return [(qubit,1),]
            
        true_matrix = true_matrix.normalized()
        true_out = matrix_to_qubit(true_matrix)

        false_matrix = false_matrix.normalized()
        false_out = matrix_to_qubit(false_matrix)
            
        return [(true_out, prob_true), (false_out, prob_false)]
    else:
        raise NotImplementedError("This function can't handle non-sympy" +\
                                  "matrix formats yet")
                                     
def measure_partial_oneshot(qubit, bit, format='sympy'):
    import random
    m = qubit_to_matrix(qubit, format)
    
    if format == 'sympy':
        m = m.normalized()
        size = max(m.shape) #max of shape to account for bra or ket
        nqubits = int(math.log(size,2)+.1)
        bit_mask = 1<<bit #when bit_mask is '&' with an array index,
        #It will determine if the bit is true in that state
        
        #break the matrix into halves one where specified bit is true, other not
        true_matrix  = zeros((2**nqubits, 1))
        false_matrix = zeros((2**nqubits, 1)) 
        for i in range(2**nqubits):
            if bit_mask&i:
                true_matrix[i] = m[i]
            else:
                false_matrix[i] = m[i]

        #calculate probability of finding the specified bit as true
        prob_true = 0
        for item in true_matrix*true_matrix.H:
            prob_true += item
        
        if random.random() < prob_true:
            return_matrix = true_matrix.normalized()
        else:
            return_matrix = false_matrix.normalized()
        
        return_matrix = matrix_to_qubit(return_matrix)

        return return_matrix        
    
def measure_all_oneshot(qubit, format='sympy'):
    import random
    m = qubit_to_matrix(qubit)
    
    if format == 'sympy':
        m = m.normalized()
        random_number = random.random()
        total = 0
        result = 0
        for i in m:
            total += i*i.conjugate()
            if total > random_number:
                break
            result += 1
        return Qubit(IntQubit(result, int(math.log(max(m.shape),2)+.1)))
    else:
        raise NotImplementedError("This function can't handle non-sympy" +\
                                  "matrix formats yet")            
