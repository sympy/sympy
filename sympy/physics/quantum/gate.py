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
* Fix Toffoli Gate to be a controlled gate.
"""

from itertools import chain

from sympy import Mul, Pow, Integer, I, pi, Matrix, Rational
from sympy.core.numbers import Number
from sympy.core.containers import Tuple
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices import matrices
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.utilities.iterables import all

from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import ComplexSpace, HilbertSpaceError
from sympy.physics.quantum.operator import UnitaryOperator
from sympy.physics.quantum.tensorproduct import matrix_tensor_product
from sympy.physics.quantum.matrixcache import matrix_cache

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

    def _represent_ZGate(self, basis, **options):
        format = options.pop('format','sympy')
        nqubits = options.pop('nqubits',0)
        if nqubits == 0:
            raise QuantumError('The number of qubits must be given as nqubits.')

        # Make sure we have enough qubits for the gate.
        if nqubits < self.min_qubits:
            raise QuantumError(
                'The number of qubits %r is too small for the gate.' % nqubits
            )

        target_matrix = self.get_target_matrix(format)
        targets = self.targets
        if isinstance(self, CGate):
            controls = self.controls
        else:
            controls = []
        m = represent_zbasis(
            controls, targets, target_matrix, nqubits, format
        )
        return m

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

    def decompose(self, **options):
        return self

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
        return matrix_cache.get_matrix('H', format)


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
        return matrix_cache.get_matrix('X', format)


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
        return matrix_cache.get_matrix('Y', format)


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
        return matrix_cache.get_matrix('Z', format)


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
        return matrix_cache.get_matrix('S', format)


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
        return matrix_cache.get_matrix('T', format)

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
        return matrix_cache.get_matrix('SWAP', format)


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

    def get_target_matrix(self, format='sympy'):
        return matrix_cache.get_matrix('X', format)


# Aliases for gate names.
CNOT = CNotGate
SWAP = SwapGate
TOFFOLI = ToffoliGate

#-----------------------------------------------------------------------------
# Represent
#-----------------------------------------------------------------------------


def _np_tensor_product(*product):
    """numpy version of tensor product of multiple arguments."""
    import numpy as np
    answer = product[0]
    for item in product[1:]:
        answer = np.kron(answer, item)
    return answer


def _np_eye(n):
    """numpy version of complex eye."""
    import numpy as np
    return np.matrix(np.eye(n, dtype='complex'))


def _sp_tensor_product(*product):
    """scipy.sparse version of tensor product of multiple arguments."""
    from scipy import sparse
    answer = product[0]
    for item in product[1:]:
        answer = sparse.kron(answer, item)
    # The final matrices will just be multiplied, so csr is a good final
    # sparse format.
    return sparse.csr_matrix(answer)


def _sp_eye(n):
    """scipy.sparse version of complex eye."""
    from scipy import sparse
    return sparse.eye(n, n, dtype='complex')


def _get_represent_utils(format='sympy'):
    """Get the version of eye and tensor_product for a given format."""
    if format == 'sympy':
        e = matrices.eye
        tp = matrix_tensor_product
    elif format == 'numpy':
        e = _np_eye
        tp = _np_tensor_product
    elif format == 'scipy.sparse':
        e = _sp_eye
        tp = _sp_tensor_product
    else:
        raise NotImplementedError('Invalid format: %r' % format) 
    return e, tp


def represent_zbasis(controls, targets, target_matrix, nqubits, format='sympy'):
    """Represent a gate with controls, targets and target_matrix.

    This function does the low-level work of representing gates as matrices
    in the standard computational basis (ZGate). Currently, we support two
    main cases:

    1. One target qubit and no control qubits.
    2. One target qubits and multiple control qubits.

    For the base of multiple controls, we use the following expression [1]:

    1_{2**n} + (|1><1|)^{(n-1)} x (target-matrix - 1_{2})

    Parameters
    ----------
    controls : list, tuple
        A sequence of control qubits.
    targets : list, tuple
        A sequence of target qubits.
    target_matrix : sympy.Matrix, numpy.matrix, scipy.sparse
        The matrix form of the transformation to be performed on the target
        qubits.  The format of this matrix must match that passed into
        the `format` argument.
    nqubits : int
        The total number of qubits used for the representation.
    format : str
        The format of the final matrix ('sympy', 'numpy', 'scipy.sparse').

    Examples
    --------

    References
    ----------
    [1] http://www.johnlapeyre.com/qinf/qinf_html/node6.html.
    """
    controls = [int(x) for x in controls]
    targets = [int(x) for x in targets]
    nqubits = int(nqubits)

    # This checks for the format as well.
    eye, tensor_product = _get_represent_utils(format)
    up = matrix_cache.get_matrix('up', format)
    eye2 = matrix_cache.get_matrix('eye2', format)

    # Plain single qubit case
    if len(controls) == 0 and len(targets) == 1:
        product = []
        bit = targets[0]
        # Fill product with [I1,Gate,I2] such that the unitaries,
        # I, cause the gate to be applied to the correct Qubit
        if bit != nqubits-1:
            product.append(eye(2**(nqubits-bit-1)))
        product.append(target_matrix)
        if bit != 0:
            product.append(eye(2**bit))
        return tensor_product(*product)

    # Single target, multiple controls.
    elif len(targets) == 1 and len(controls) >= 1:
        target =  targets[0]

        # Build the non-trivial part.
        product2 = []
        for i in range(nqubits):
            product2.append(eye(2))
        for control in controls:
            product2[nqubits-1-control] = up
        product2[nqubits-1-target] = target_matrix - eye2

        return eye(2**nqubits) + tensor_product(*product2)

    # Multi-target, multi-control is not yet implemented.
    else:
        raise NotImplementedError(
            'The representation of multi-target, multi-control gates '
            'is not implemented.'
        )


#-----------------------------------------------------------------------------
# Gate manipulation functions.
#-----------------------------------------------------------------------------


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
                    newargs = newargs + (PhaseGate(circuit.args[i].base.args[0][0])**\
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
                    #meaning they do not affect ANY of the same targets
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


#-----------------------------------------------------------------------------
# Utility functions
#-----------------------------------------------------------------------------


def zx_basis_transform(self, format='sympy'):
    """Transformation matrix from Z to X basis."""
    return matrix_cache.get_matrix('ZX', format)


def zy_basis_transform(self, format='sympy'):
    """Transformation matrix from Z to Y basis."""
    return matrix_cache.get_matrix('ZY', format)

