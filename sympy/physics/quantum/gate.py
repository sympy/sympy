"""An implementation of qubits and gates acting on them.

Todo:
* Optimize Gate._apply_operators_Qubit to remove the creation of many 
  intermediate Qubit objects.
* Optimize the get_target_matrix by using slots to precompute the matrices.
* Get UGate to work with either sympy/numpy matrices and output either
  format. This should also use the matrix slots.
* Add commutation relationships to all operators and use this in gate_sort.
* Get represent working and test.
* Test apply_operators.
"""

from itertools import chain

from sympy import Mul, Pow, Integer, I, pi, Matrix, Rational
from sympy.core.numbers import Number
from sympy.core.containers import Tuple
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.utilities.iterables import all

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import ComplexSpace, HilbertSpaceError
from sympy.physics.quantum.operator import UnitaryOperator
from sympy.physics.quantum.tensorproduct import matrix_tensor_product

__all__ = [
    'Gate',
    'CGate',
    'UGate',
    'OneQubitGate',
    'TwoQubitGate',
    'HadamardGate',
    'XGate',
    'YGate',
    'ZGate',
    'TGate',
    'PhaseGate',
    'SwapGate',
    'CNotGate',
    'ToffoliGate',
    # Aliased gate names
    'U',
    'CNOT',
    'SWAP',
    'TOFFOLI',
    'H',
    'X',
    'Y',
    'Z',
    'T',
    'S',
    'Phase'
]

sqrt2_inv = Pow(2, Rational(-1,2), evaluate=False)

#-----------------------------------------------------------------------------
# Gate Super-Classes
#-----------------------------------------------------------------------------

def _validate_targets_controls(tandc):
    tandc = list(tandc)
    # Check for integers
    for bit in tandc:
        if not bit.is_Integer:
            raise TypeError('Integer expected, got: %r' % tandc[bit])
    # Detect duplicates
    if len(list(set(tandc))) != len(tandc):
        raise QuantumError(
            'Target/control qubits in a gate cannot be duplicated'
        )

class Gate(UnitaryOperator):
    """Non-controlled unitary gate operator that acts on qubits.

    This is a general abstract gate that needs to be subclassed to do anything
    useful.

    Parameters
    ----------
    label : tuple, int
        A list of the target qubits (as ints) that the gate will apply to.

    Examples
    --------


    """

    _label_separator = ','

    gate_name = u'G'
    gate_name_latex = u'G'

    #-------------------------------------------------------------------------
    # Initialization/creation
    #-------------------------------------------------------------------------

    @classmethod
    def _eval_label(cls, label):
        label = UnitaryOperator._eval_label(label) # Make label a Tuple and sympify.
        _validate_targets_controls(label)
        return label

    @classmethod
    def _eval_hilbert_space(cls, label):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(label)+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def nqubits(self):
        """The total number of qubits this gate acts on.

        For controlled gate subclasses this includes both target and control
        qubits, so that, for examples the CNOT gate acts on 2 qubits.
        """
        return len(self.targets)

    @property
    def min_qubits(self):
        """The minimum number of qubits this gate needs to act on."""
        return max(self.targets)+1

    @property
    def targets(self):
        """A tuple of target qubits."""
        return self.label

    #-------------------------------------------------------------------------
    # Gate methods
    #-------------------------------------------------------------------------

    def get_target_matrix(self, format='sympy'):
        """The matrix rep. of the target part of the gate.

        Parameters
        ----------
        format : str
            The format string ('sympy','numpy', etc.)
        """
        raise NotImplementedError('get_target_matrix is not implemented in Gate.')

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    def _apply_operator_Qubit(self, qubits, **options):
        """Apply this gate to a Qubit."""

        # Check number of qubits this gate acts on.
        if qubits.nqubits < self.min_qubits:
            raise QuantumError(
                'Gate needs a minimum of %r qubits to act on, got: %r' %\
                    (self.min_qubits, qubits.nqubits)
            )

        # If the controls are not met, just return
        if isinstance(self, CGate):
            if not self.eval_controls(qubits):
                return qubits

        targets = self.targets
        target_matrix = self.get_target_matrix(format='sympy')

        # Find which column of the target matrix this applies to.
        column_index = 0
        n = 1
        for target in targets:
            column_index += n*qubits[target]
            n = n<<1
        column = target_matrix[:,int(column_index)]
    
        # Now apply each column element to the qubit.
        result = 0
        for index in range(column.rows):
            # TODO: This can be optimized to reduce the number of Qubit
            # creations. We should simply manipulate the raw list of qubit
            # values and then build the new Qubit object once.
            # Make a copy of the incoming qubits.
            new_qubit = qubits.__class__(*qubits.args)
            # Flip the bits that need to be flipped.
            for bit in range(len(targets)):
                if new_qubit[targets[bit]] != (index>>bit)&1:
                    new_qubit = new_qubit.flip(targets[bit])
            # The value in that row and column times the flipped-bit qubit
            # is the result for that part.
            result += column[index]*new_qubit
        return result

    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    # def _represent(self, basis, format='sympy'):
    #     if isinstance(basis, Gate):
    #         basis_size = 1
    #     elif isinstance(basis, Pow) and isinstance(basis.base, Gate):
    #         basis_size = basis.exp
    #     else:
    #         raise QuantumError(
    #             'Basis must be a gate operator, or tensor product of gate operators.'
    #         )
    # 
    #     if basis_size < self.min_qubits:
    #         raise QuantumError(
    #             'The basis given is too small for the given Gate objects.'
    #         )
    # 
    #     gate = self.matrix
    #     if isinstance(basis, Gate):
    #         return gate
    #     else:
    #         m = represent_hilbert_space(
    #             gate, basis_size, self.label, format
    #         )
    #         return m

    #-------------------------------------------------------------------------
    # Print methods
    #-------------------------------------------------------------------------

    def _print_contents(self, printer, *args):
        label = self._print_label(printer, *args)
        return '%s(%s)' % (self.gate_name, label)

    def _print_contents_pretty(self, printer, *args):
        a = stringPict(unicode(self.gate_name))
        b = self._print_label_pretty(printer, *args)
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _latex(self, printer, *args):
        label = self._print_label(printer, *args)
        return '%s_{%s}' % (self.gate_name_latex, label)


class CGate(Gate):
    """A general unitary gate with control qubits.

    A general control gate applies a target gate to a set of targets if all
    of the control qubits have a particular values (set by
    ``CGate.control_value``).

    Parameters
    ----------
    label : tuple
        The label in this case has the form (controls, gate), where controls
        is a tuple/list of control qubits (as ints) and gate is a ``Gate``
        instance that is the target operator.

    Examples
    --------

    """

    gate_name = u'C'
    gate_name_latex = u'C'

    # The values this class controls for.
    control_value = Integer(1)

    #-------------------------------------------------------------------------
    # Initialization
    #-------------------------------------------------------------------------

    # CGate(((0,1),Gate(0)))

    @classmethod
    def _eval_label(cls, label):
        # _eval_label has the right logic for the controls argument.
        controls = UnitaryOperator._eval_label(label[0])
        gate = label[1]
        _validate_targets_controls(chain(controls,gate.targets))
        return Tuple(controls, gate)

    @classmethod
    def _eval_hilbert_space(cls, label):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**max(max(label[0])+1,label[1].min_qubits)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def nqubits(self):
        """The total number of qubits this gate acts on.

        For controlled gate subclasses this includes both target and control
        qubits, so that, for examples the CNOT gate acts on 2 qubits.
        """
        return len(self.targets)+len(self.controls)

    @property
    def min_qubits(self):
        """The minimum number of qubits this gate needs to act on."""
        return max(max(self.controls),max(self.targets))+1

    @property
    def targets(self):
        """A tuple of target qubits."""
        return self.gate.targets

    @property
    def controls(self):
        """A tuple of control qubits."""
        return self.label[0]

    @property
    def gate(self):
        """The non-controlled gate that will be applied to the targets."""
        return self.label[1]

    #-------------------------------------------------------------------------
    # Gate methods
    #-------------------------------------------------------------------------

    def get_target_matrix(self, format='sympy'):
        return self.gate.get_target_matrix(format)

    def eval_controls(self, qubit):
        """Return True/False to indicate if the controls are satisfied."""
        return all([qubit[bit]==self.control_value for bit in self.controls])

    #-------------------------------------------------------------------------
    # Print methods
    #-------------------------------------------------------------------------

    def _print_controls(self, printer, *args):
        result = []
        for item in self.controls:
            result.append(printer._print(item, *args))
        return self._label_separator.join(result)

    def _print_controls_pretty(self, printer, *args):
        pform = printer._print(self.controls[0], *args)
        for item in self.controls[1:]:
            pform = prettyForm(*pform.right((self._label_separator)))
            nextpform = printer._print(item, *args)
            pform = prettyForm(*pform.right((nextpform)))
        return pform

    def _print_contents(self, printer, *args):
        controls = self._print_controls(printer, *args)
        gate = printer._print(self.gate, *args)
        return '%s(((%s),%s))' %\
            (self.gate_name, controls, gate)

    def _print_contents_pretty(self, printer, *args):
        controls = self._print_controls_pretty(printer, *args)
        gate = printer._print(self.gate)
        gate_name = stringPict(unicode(self.gate_name))
        top = stringPict(*controls.left(' '*gate_name.width()))
        bot = stringPict(*gate_name.right(' '*controls.width()))
        first = prettyForm(binding=prettyForm.POW, *bot.below(top))
        gate = prettyForm(*gate.parens(left='(', right=')'))
        final = prettyForm(*first.right((gate)))
        return final

    def _latex(self, printer, *args):
        controls = self._print_controls(printer, *args)
        gate = printer._print(self.gate, *args)
        return r'%s_{%s}{\left(%s\right)}' %\
            (self.gate_name_latex, controls, gate)


class UGate(Gate):
    """General gate specified by a set of targets and a target matrix.

    Parameters
    ----------
    label : tuple
        A tuple of the form (targets, U), where targets is a tuple of the
        target qubits and U is a unitary matrix with dimension of
        len(targets).
    """
    gate_name = u'U'
    gate_name_latex = u'U'

    #-------------------------------------------------------------------------
    # Initialization
    #-------------------------------------------------------------------------

    @classmethod
    def _eval_label(cls, label):
        # Gate._eval_label has the right logic.
        targets = Gate._eval_label(label[0])
        mat = label[1]
        if not isinstance(mat, Matrix):
            raise TypeError('Matrix expected, got: %r' % mat)
        dim = 2**len(targets)
        if not all([dim == shape for shape in mat.shape]):
        # if (dim != mat.shape[0]) or (dim != mat.shape[1]):
            raise IndexError(
                'Number of targets must match the matrix size: %r %r' %\
                (targets, mat)
            )
        return Tuple(targets, mat)

    @classmethod
    def _eval_hilbert_space(cls, label):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(label[0])+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def targets(self):
        """A tuple of target qubits."""
        return self.label[0]

    #-------------------------------------------------------------------------
    # Gate methods
    #-------------------------------------------------------------------------

    def get_target_matrix(self, format='sympy'):
        """The matrix rep. of the target part of the gate.

        Parameters
        ----------
        format : str
            The format string ('sympy','numpy', etc.)
        """
        return self.label[1]

    #-------------------------------------------------------------------------
    # Print methods
    #-------------------------------------------------------------------------

    def _print_targets(self, printer, *args):
        result = []
        for item in self.targets:
            result.append(printer._print(item, *args))
        return self._label_separator.join(result)

    def _print_targets_pretty(self, printer, *args):
        pform = printer._print(self.targets[0], *args)
        for item in self.targets[1:]:
            pform = prettyForm(*pform.right((self._label_separator)))
            nextpform = printer._print(item, *args)
            pform = prettyForm(*pform.right((nextpform)))
        return pform

    def _print_contents(self, printer, *args):
        targets = self._print_targets(printer, *args)
        return '%s(%s)' % (self.gate_name, targets)

    def _print_contents_pretty(self, printer, *args):
        targets = self._print_targets_pretty(printer, *args)
        gate_name = stringPict(unicode(self.gate_name))
        top = stringPict(*targets.left(' '*gate_name.width()))
        bot = stringPict(*gate_name.right(' '*targets.width()))
        result = prettyForm(binding=prettyForm.POW, *bot.below(top))
        return result

    def _latex(self, printer, *args):
        targets = self._print_targets(printer, *args)
        return r'%s_{%s}' % (self.gate_name_latex, targets)


class OneQubitGate(Gate):
    """A single qubit unitary gate base class."""

    nqubits = Integer(1)


class TwoQubitGate(Gate):
    """A two qubit unitary gate base class."""

    nqubits = Integer(2)


# Aliases for gate names.
U = UGate

#-----------------------------------------------------------------------------
# Single Qubit Gates
#-----------------------------------------------------------------------------


class HadamardGate(OneQubitGate):
    """The single qubit Hadamard gate.

    Parameters
    ----------
    target : int
        The target qubit this gate will apply to.

    Examples
    --------

    """
    gate_name = u'H'
    gate_name_latex = u'H'
        
    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return sqrt2_inv*Matrix([[1, 1], [1, -1]])
        raise NotImplementedError('Invalid format: %r' % format)


class XGate(OneQubitGate):
    """The single qubit X, or NOT, gate.

    Parameters
    ----------
    target : int
        The target qubit this gate will apply to.

    Examples
    --------

    """
    gate_name = u'X'
    gate_name_latex = u'X'

    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[0, 1], [1, 0]])
        raise NotImplementedError('Invalid format: %r' % format)


class YGate(OneQubitGate):
    """The single qubit Y gate.

    Parameters
    ----------
    target : int
        The target qubit this gate will apply to.

    Examples
    --------

    """
    gate_name = u'Y'
    gate_name_latex = u'Y'
        
    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[0, complex(0,-1)], [complex(0,1), 0]])
        raise NotImplementedError('Invalid format: %r' % format)


class ZGate(OneQubitGate):
    """The single qubit Z gate.

    Parameters
    ----------
    target : int
        The target qubit this gate will apply to.

    Examples
    --------

    """
    gate_name = u'Z'
    gate_name_latex = u'Z'

    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[1, 0], [0, -1]])
        raise NotImplementedError('Invalid format: %r' % format)


class PhaseGate(OneQubitGate):
    """The single qubit phase gate.

    This gate rotates the phase of the state by pi/2 if the state is |1> and
    does nothing if the state is |0>.

    Parameters
    ----------
    target : int
        The target qubit this gate will apply to.

    Examples
    --------

    """
    gate_name = u'S'
    gate_name_latex = u'S'

    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[1, 0], [0, complex(0,1)]])
        raise NotImplementedError('Invalid format: %r' % format)


class TGate(OneQubitGate):
    """The single qubit pi/8 gate.

    This gate rotates the phase of the state by pi/4 if the state is |1> and
    does nothing if the state is |0>.

    Parameters
    ----------
    target : int
        The target qubit this gate will apply to.

    Examples
    --------

    """
    gate_name = u'T'
    gate_name_latex = u'T'

    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[1, 0], [0, exp(I*pi/4)]])
        raise NotImplementedError('Invalid format: %r' % format)

# Aliases for gate names.
H = HadamardGate
X = XGate
Y = YGate
Z = ZGate
T = TGate
Phase = S = PhaseGate


#-----------------------------------------------------------------------------
# 2 Qubit Gates
#-----------------------------------------------------------------------------


class CNotGate(CGate, TwoQubitGate):
    """Two qubit controlled-NOT.

    This gate performs the NOT or X gate on the target qubit if the control
    qubits all have the value 1.

    Parameters
    ----------
    label : tuple
        A tuple of the form (control, target).

    Examples
    --------

    """
    gate_name = 'CNOT'
    gate_name_latex = u'CNOT'

    #-------------------------------------------------------------------------
    # Initialization
    #-------------------------------------------------------------------------

    @classmethod
    def _eval_label(cls, label):
        label = Gate._eval_label(label)
        return label

    @classmethod
    def _eval_hilbert_space(cls, label):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(label)+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def min_qubits(self):
        """The minimum number of qubits this gate needs to act on."""
        return max(self.label)+1

    @property
    def targets(self):
        """A tuple of target qubits."""
        return Tuple(self.label[1])

    @property
    def controls(self):
        """A tuple of control qubits."""
        return Tuple(self.label[0])

    @property
    def gate(self):
        """The non-controlled gate that will be applied to the targets."""
        return XGate(self.label[1])

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    # The default printing of Gate works better than those of CGate, so we
    # go around the overridden methods in CGate.

    def _print_contents(self, printer, *args):
        return Gate._print_contents(self, printer, *args)

    def _print_contents_pretty(self, printer, *args):
        return Gate._print_contents_pretty(self, printer, *args)

    def _latex(self, printer, *args):
        return Gate._latex(self, printer, *args)


class SwapGate(TwoQubitGate):
    """Two qubit SWAP gate.

    This gate swap the values of the two qubits.

    Parameters
    ----------
    label : tuple
        A tuple of the form (target1, target2).

    Examples
    --------

    """
    gate_name = 'SWAP'
    gate_name_latex = u'SWAP'

    def get_target_matrix(self, format='sympy'):
        if format == 'sympy':
            return Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
        raise NotImplementedError('Invalid format: %r' % format)


class ToffoliGate(Gate):
    """The ToffoliGate, also known as the double-controlled-NOT gate.

    Flips the third Qubit (the target) contingent on the first and
    second qubits (the controls) being 1.

    Parameters
    ----------
    label : tuple
        A tuple of the form (control1, control2, target).

    Examples
    --------

    """

    nqubits = Integer(3)

    gate_name = u'TOFFOLI'
    gate_name_latex = u'TOFFOLI'

    #-------------------------------------------------------------------------
    # Initialization
    #-------------------------------------------------------------------------

    # Toffoli((0,1,2))

    @classmethod
    def _eval_label(cls, label):
        label = Gate._eval_label(label)
        return label

    @classmethod
    def _eval_hilbert_space(cls, label):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(label)+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def min_qubits(self):
        """The minimum number of qubits this gate needs to act on."""
        return max(self.label)+1

    @property
    def targets(self):
        """A tuple of target qubits."""
        return Tuple(self.label[2])

    @property
    def controls(self):
        """A tuple of control qubits."""
        return self.label[:2]

    @property
    def gate(self):
        """The non-controlled gate that will be applied to the targets."""
        return XGate(self.label[2])

# Aliases for gate names.
CNOT = CNotGate
SWAP = SwapGate
TOFFOLI = ToffoliGate

#-----------------------------------------------------------------------------
# Utility functions
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


def zx_basis_transform(self, format='sympy'):
    """Transformation matrix from Z to X basis."""
    if format == 'sympy':
        return sqrt2_inv*Matrix([[1,1],[1,-1]])
    raise NotImplementedError('Invalid format: %r' % format)


def zy_basis_transform(self, format='sympy'):
    """Transformation matrix from Z to Y basis."""
    if format == 'sympy':
        return Matrix([[I,0],[0,-I]])
    raise NotImplementedError('Invalid format: %r' % format)
