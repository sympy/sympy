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
                                  
def measure_partial(qubit, bits, format='sympy'):
    """Does a partial ensemble measure on the specifed qubit 
       Qubit is the state of the system
       bits is an array or tuple of bits to measure
       
       Examples
       ==========
       state = Qubit(0,1)+Qubit(1,0)
       measure_partial(state,(0,))
       
    """
    m = qubit_to_matrix(qubit, format)

    if format == 'sympy':
        result = []
        m = m.normalized()
        possible_outcomes = __get_possible_outcomes(m, bits)

        #form output from function
        output = []
        for outcome in possible_outcomes:
            #calculate probability of finding the specified bits with given values
            prob_of_outcome = 0
            prob_of_outcome += (outcome.H*outcome)[0]
            
            #If the output has a chance, append it to output with found probability    
            if prob_of_outcome != 0:    
                output.append((matrix_to_qubit(outcome.normalized()),\
                prob_of_outcome))
            
        return output
    else:
        raise NotImplementedError("This function can't handle non-sympy" +\
                                  "matrix formats yet")
                                     
def measure_partial_oneshot(qubit, bits, format='sympy'):
    import random
    m = qubit_to_matrix(qubit, format)
    
    if format == 'sympy':
        result = []
        m = m.normalized()
        possible_outcomes = __get_possible_outcomes(m, bits) 

        #form output from function
        output = []
        random_number = random.random()
        total_prob = 0
        for outcome in possible_outcomes:
            #calculate probability of finding the specified bits with given values
            total_prob += (outcome.H*outcome)[0]
            if total_prob >= random_number:
                return matrix_to_qubit(outcome.normalized())
    else:
        raise NotImplementedError("This function can't handle non-sympy" +\
                                  "matrix formats yet")
def __get_possible_outcomes(m, bits):
    """
        Function inputs:
            m: the matrix representing the state of the system
            bits: a tuple or list with the bits that will be measured
        
        Function outputs:
            The list of possible states which can occur given this measurement.
            These are un-normalized so we can derive the probability of finding
            this state by taking the inner product with itself
             
        This is filled with loads of dirty binary tricks...You have been warned
    
    """
    size = max(m.shape) #max of shape to account for bra or ket
    nqubits = int(math.log(size,2)+.1) #number of qubits possible

    #Make the output states and put in output_matrices, nothing in them now.
    #Each state will represent a possible outcome of the measurement
    #Thus, output_matrices[0] is the matrix which we get when all measured bits
    #return 0. and output_matrices[1] is the matrix for only the 0th bit being true 
    output_matrices = []
    for i in range(1<<len(bits)):
        output_matrices.append(zeros((2**nqubits, 1)))    

    
    #Bitmasks will help sort how to determine possible outcomes.
    #When the bit mask is and-ed with a matrix-index, It will determine
    #it will determine which state that index belongs to
    bit_masks = [] 
    for bit in bits:
        bit_masks.append(1<<bit)
     


    #make possible outcome states
    for i in range(2**nqubits):
        trueness = 0 #This tells us to which output_matrix this value belongs
        #Find trueness
        for j in range(len(bit_masks)):
            if i&bit_masks[j]:
                trueness += j+1
        #put the value in the correct output matrix        
        output_matrices[trueness][i] = m[i]    
    return output_matrices
    
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
