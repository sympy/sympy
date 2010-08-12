from sympy.physics.hilbert import l2, HilbertSpaceException
from sympy.physics.quantum import BasisSet, Operator, Representable, represent, OuterProduct, Ket
from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify, Matrix, elementary
from sympy.core.numbers import *
from sympy.core.basic import S, sympify
from sympy.core.function import Function
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices.matrices import Matrix, eye
from sympy.core.symbol import Symbol, symbols
from sympy.physics.qmul import QMul
from sympy.physics.qadd import QAdd
from sympy.physics.qpow import QPow
from sympy.physics.quantumbasic import QuantumError
import math

#-----------------------------------------------------------------------------
# BasisSets for Qbit
#-----------------------------------------------------------------------------
class QbitZBasisSet(BasisSet):
    """
        ZBasisSet formed using the eigenvectors of the Pauli-Z Matrix
    """
    __slots__ = ['hilbert_space', 'qbit_number']

    def __new__(cls, nqbits):
        obj = BasisSet.__new__(cls, 2**nqbits)
        obj.hilbert_space = l2(2)**nqbits
        obj.qbit_number = nqbits
        return obj

class QbitXBasisSet(BasisSet):
    """
        XBasisSet formed using the eigenvectors of the Pauli-X Matrix
    """
    __slots__ = ['hilbert_space', 'qbit_number']

    def __new__(cls, nqbits):
        obj = BasisSet.__new__(cls, 2**nqbits)
        obj.hilbert_space = l2(2)**nqbits
        obj.qbit_number = nqbits
        return obj

class QbitYBasisSet(BasisSet):
    """
        XBasisSet formed using the eigenvectors of the Pauli-X Matrix
    """
    __slots__ = ['hilbert_space', 'qbit_number']

    def __new__(cls, nqbits):
        obj = BasisSet.__new__(cls, 2**nqbits)
        obj.hilbert_space = l2(2)**nqbits
        obj.qbit_number = nqbits
        return obj

#-----------------------------------------------------------------------------
# Qbit Classes
#-----------------------------------------------------------------------------

class Qbit(Ket):
    """
        Represents a single definite eigenstate of the ZBasisSet
        Object can be instantiated by:
            - *args with each representing the value of a single Qbit
            - A single decimal value. Code will convert to binary using least number of bits
            - A decimal value and the number of bits you wish to express it in (The size of the Hilbert Space)

        Qbit object can have individual bits flipped such that a one becomes a zero and vice versa

        Qbit object is also callable allowing user to pick out the value of a particular bit

        Can also set outDecimal class variable which causes the state to present itself in decimal

        All of this is seen in the doc test:

        >>> from sympy.physics.qbit import Qbit
        >>> Qbit(0,0,0)
        |'000'>
        >>> Qbit(5)
        |'101'>
        >>> a = Qbit(5,4)
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
        >>> Qbit.outDecimal = True
        >>> a
        |'5'>
        >>> a.to_string()
        '5'
    """
    outDecimal = False
    def __new__(cls, *args, **options):
        import math
        #If they just give us one number, express it in the least number of bits possible
        if args[0] > 1 and len(args) == 1:
            array = [(args[0]>>i)&1 for i in reversed(range(int(math.ceil(math.log(args[0], 2)+.01)+.001)))]
            array = sympify(array)
            obj = Expr.__new__(cls, *array)
            obj.evaluates = cls
            obj.hilbert_space = l2(2)
            return obj
        #if they give us two numbers, the second number is the number of bits on which it is expressed)
        #Thus, Qbit(0,5) == |00000>. second argument can't be one becuase of intersection and uslessesness of possibility
        elif len(args) == 2 and args[1] > 1:
            array = [(args[0]>>i)&1 for i in reversed(range(args[1]))]
            array = sympify(array)
            obj = Expr.__new__(cls, *array)
            obj.evaluates = cls
            obj.hilbert_space = l2(2)
            return obj
        for element in args:
            if not (element == 1 or element == 0):
                raise QuantumError("Values must be either one or zero")
        args = sympify(args)
        obj = Expr.__new__(cls, *args)
        obj.evaluates = cls
        obj.hilbert_space = l2(2)
        return obj

    @property
    def dimension(self):
        """
            returns the number of qbits in object
        """
        return len(self.args)

    def __len__(self):
        return self.dimension

    def __getitem__(self, bit):
        if bit > self.dimension - 1:
            raise QuantumError()
        return self.args[int(self.dimension-bit-1)]

    def _print_name(self, printer, *args):
        string = self.to_string()
        return printer._print(string, *args)

    def _print_name_pretty(self, printer, *args):
        string = self.to_string()
        pform = printer._print(string, *args)
        return pform

    def to_string(self):
        """
            Returns a string representing the *args tuple
        """
        if Qbit.outDecimal:
            number = 0
            n = 1
            for i in reversed(self.args):
                number += n*i
                n = n<<1
            string = str(number)
        else:
            string = ""
            for i in self.args:
                string = string + str(i)
        return string

    def _sympyrepr(self, printer, *args):
        return "%s%s" % (printer._print(self.__class__.__name__, *args), printer._print(str(self.args[0])))

    def flip(self, *args):
        """
            flips the bit(s) given in *args
        """
        newargs = list(self.args)
        for i in args:
            bit = int(self.dimension-i-1)
            if newargs[bit] == 1:
                newargs[bit] = 0
            else:
                newargs[bit] = 1
        return Qbit(*newargs)

    def _represent_QbitZBasisSet(self, basis, **options):
        #make sure HilbertSizes of basis and qbit match
        if basis.hilbert_space != l2(2)**self.dimension:
            raise HilbertSpaceException("Basis and Qbit dimensions do not match!")
        n = 1
        definiteState = 0
        args = self.args
        for it in reversed(args):
            definiteState += n*it
            n = n*2
        result = [0 for x in range(2**(self.dimension))]
        result[int(definiteState)] = 1
        return Matrix(result)

    @property
    def XBasisTransform(self):
        #Transform matrix from ZBasis to XBasis and back again
        return 1/sqrt(2)*Matrix([[1,1],[1,-1]])

    @property
    def YBasisTransform(self):
        #Transform matrix from ZBasis to YBasis and back again
        return Matrix([[ImaginaryUnit(),0],[0,-ImaginaryUnit()]])

class QbitX(Qbit):
    def __new__(cls, *args):
        for element in args:
            if not (element == '+' or element == '-'):
                raise QuantumError("Values must be either + or -")
        return Expr.__new__(cls, *args, **{'commutative': False})

#-----------------------------------------------------------------------------
# Gate Super-Classes
#-----------------------------------------------------------------------------
class Gate(Operator):
    """
        The superclass of gate operators that acts on qubit(s)
            - *args stores the information about which qbits it will act on
            - minimumdimension is the least size a set fo qbits must be for the gate to be applied
            - input number is the number of qbits the gate acts on (e.g. a three qbit gate acts on three qbits
            - name returns the string of the qbit
            - apply applies gate onto qbits
            - _represent_Qbit?BasisSet represents the gates matrix in the given basis
    """

    def __new__(cls, *args):
        args = sympify(args)
        obj = Expr.__new__(cls, *args)
        if obj.inputnumber != len(args):
            num = obj.inputnumber
            raise QuantumError("This gate applies to %d qbits" % (num))
        for i in range(len(args)):
            if args[i] in (args[:i] + args[i+1:]):
                raise QuantumError("Can't have duplicate control and target bits!")
        return obj

    @property
    def hilbert_space(self):
        return l2(2)

    @property
    def matrix(self):
        raise NotImplementedError("matrixRep Not implemented")

    @property
    def minimumdimension(self):
        return max(self.args)

    @property
    def inputnumber(self):
        #number of inputs is determined by the size of the matrix
        import math
        mat = self.matrix
        return int(math.log((mat).cols,2)+.5)

    def _apply(self, qbits, mat, args):
        assert isinstance(qbits, Qbit), "can only apply self to qbits"
        # check number of qbits this gate acts on
        if self.minimumdimension >= qbits.dimension:
            raise HilbertSpaceException()
        if isinstance(qbits, Qbit):
            #find which column of the matrix this qbit applies to
            column_index = 0
            n = 1
            for element in args:
                column_index += n*qbits[element]
                n = n<<1
            column = mat[:,int(column_index)]
            #now apply each column element to qbit
            result = 0
            for index in range(len(column.tolist())):
                new_qbit = Qbit(*qbits.args)
                #flip the bits that need to be flipped
                for bit in range(len(args)):
                    if new_qbit[args[bit]] != (index>>bit)&1:
                        new_qbit = new_qbit.flip(args[bit])
                #the value in that row and column times the flipped-bit qbit is the result for that part
                result += column[index]*new_qbit
        else:
            raise QuantumError("can't apply to object that is not a qbit")
        return result

    def _apply_QbitZBasisSet(self, qbits):
        #switch qbit basis and matrix basis when fully implemented
        mat = self.matrix
        args = [self.args[i] for i in reversed(range(len(self.args)))]
        return self._apply(qbits, mat, args)

    def _sympyrepr(self, printer, *args):
        return "%s(%s)" %  (printer._print(self.__class__.__name__, *args), printer._print(self.args, *args))

    def _represent_QbitZBasisSet(self, basis, format = 'sympy'):
        if self.minimumdimension >= basis.qbit_number:
            raise HilbertSpaceException()
        gate = self.matrix
        if basis.qbit_number  == 1:
            return gate
        else:
            m = representHilbertSpace(gate, basis.qbit_number, self.args, format)
            return m

    @property
    def name(self):
        string = "("
        arguments = len(self.args)
        for i in range(arguments):
            string = string + str(self.args[i])
            if i != arguments - 1:
                string = string + ","
        string = string + ")"
        return self.__class__.__name__ + string

    @property
    def XBasisTransform(self):
        #Transform matrix from ZBasis to XBasis and back again
        return 1/sqrt(2)*Matrix([[1,1],[1,-1]])

    @property
    def YBasisTransform(self):
        #Transform matrix from ZBasis to YBasis and back again
        return Matrix([[ImaginaryUnit(),0],[0,-ImaginaryUnit()]])

class NondistributiveGate(Gate):
    """
        Superclass for all non-distributive gates (e.g. measurement)
        Gates that do not distribute should subclass off this and implement an appropriate 'measure' function
        TODO: Get Nondistributive Gates working with represent and not just apply
    """
    def __new__(cls, *args):
        args = sympify(args)
        if len(args) != 1:
            raise QuantumError("Can only measure one at a time for now")
        return Expr.__new__(cls, *args, **{'commutative': False})

    def measure(states):
        #measure is the name of non-commutative apply
        raise NotImplementedError("This type of measure has not been implemented")

#-----------------------------------------------------------------------------
# Single Qbit Gates
#-----------------------------------------------------------------------------

class HadamardGate(Gate):
    """
        An object representing a Hadamard Gate
        >>> from sympy.physics.qbit import HadamardGate, Qbit, apply_gates, QbitZBasisSet
        >>> HadamardGate(0)*Qbit(0,0)
        HadamardGate(0)*|00>
        >>> apply_gates(_)
        2**(1/2)/2*|00> + 2**(1/2)/2*|01>
        >>> from sympy.physics.quantum import represent
        >>> represent(HadamardGate(0)*Qbit(0,0), QbitZBasisSet(2))
        [2**(1/2)/2]
        [2**(1/2)/2]
        [         0]
        [         0]

    """

    @property
    def matrix(self):
        from sympy.functions.elementary.miscellaneous import sqrt
        return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))

class XGate(Gate):
    """
        An object representing a Pauli-X gate
        Also Known as a NOT Gate because it flips the value of a bit to its opposit like the classical Not Gate

        >>> from sympy.physics.qbit import Qbit, XGate, apply_gates, QbitZBasisSet
        >>> XGate(0)*Qbit(0,0)
        XGate(0)*|00>
        >>> apply_gates(_)
        |01>
        >>> from sympy.physics.quantum import represent
        >>>represent(XGate(0)*Qbit(0,0), QBitZBasisSet(2))
        >>> represent(XGate(0)*Qbit(0,0), QbitZBasisSet(2))
        [1]
        [1]
        [0]
        [0]
    """
    @property
    def matrix(self):
        return Matrix([[0, 1], [1, 0]])

class YGate(Gate):
    """
        An object representing a Pauli-Y gate
        >>> from sympy.physics.qbit import Qbit, YGate, apply_gates, QbitZBasisSet
        >>> YGate(0)*Qbit(0,1)
        YGate(0)*|01>
        >>> apply_gates(_)
        -1.0*I*|00>
        >>> from sympy.physics.quantum import represent
        >>> represent(YGate(0)*Qbit(0,0), QbitZBasisSet(2))
        [0]
        [I]
        [0]
        [0]
    """

    @property
    def matrix(self):
        return Matrix([[0, complex(0,-1)], [complex(0,1), 0]])

class ZGate(Gate):
    """
        An object representing a Pauli-Z gate:
        >>> from sympy.physics.qbit import Qbit, ZGate, apply_gates, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> ZGate(0)*Qbit(0,1)
        ZGate(0)*|01>
        >>> apply_gates(_)
        -1*|01>
        >>> represent(ZGate(0)*Qbit(0,0), QbitZBasisSet(2))
        [1]
        [0]
        [0]
        [0]

    """

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, -1]])

class PhaseGate(Gate):
    """
        An object representing a phase gate

        >>> from sympy.physics.qbit import Qbit, PhaseGate, apply_gates, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(PhaseGate(0)*Qbit(0,1), QbitZBasisSet(2))
        [0]
        [I]
        [0]
        [0]
        >>> apply_gates(PhaseGate(0)*Qbit(0,1))
        I*|01>
        >>> PhaseGate(0)*Qbit(0,1)
        PhaseGate(0)*|01>
    """

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, complex(0,1)]])

class TGate(Gate):
    """
        An object representing a pi/8 gate
        >>> from sympy.physics.qbit import Qbit, TGate, apply_gates, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> TGate(0)*Qbit(0,1)
        TGate(0)*|01>
        >>> apply_gates(_)
        exp(pi*I/4)*|01>
        >>> represent(TGate(0)*Qbit(0,1), QbitZBasisSet(2))
        [          0]
        [exp(pi*I/4)]
        [          0]
        [          0]
    """

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, exp(I*Pi()/4)]])
#-----------------------------------------------------------------------------
# 2 Qbit Gates
#-----------------------------------------------------------------------------
class RkGate(Gate):
    """
        A Controlled phase gate.
        If qbits specified in self.args[0] and self.args[1] are 1, then changes the phase of the state by e**(2*i*pi/2**k)
        k is set by the thrird argument in the input
        >>> from sympy.physics.qbit import Qbit, RkGate, apply_gates, QbitZBasisSet
        >>> RkGate(1,0,2)
        R2(1, 0)
        >>> from sympy.physics.quantum import represent
        >>> represent(_, QbitZBasisSet(2))
        [1, 0, 0, 0]
        [0, 1, 0, 0]
        [0, 0, 1, 0]
        [0, 0, 0, I]
        >>> RkGate(1,0,3)*Qbit(1,1)
        R3(1, 0)*|11>
        >>> apply_gates(_)
        exp(pi*I/4)*|11>
    """
    __slots__ = ['k']

    def __new__(cls, *args):
        obj = Expr.__new__(cls, *args[:-1], **{'commutative': False})
        if obj.inputnumber != len(args):
            num = obj.inputnumber
            raise QuantumError("This gate applies to %d qbits" % (num))
        obj.k = args[-1]
        return obj

    def _apply_QbitZBasisSet(self, qbits):
        #switch qbit basis and matrix basis when fully implemented
        mat = self.matrix
        args = [self.args[i] for i in reversed(range(2))]
        return self._apply(qbits, mat, args)

    def _sympystr(self, printer, *args):
        return "R%s(%s, %s)" % (printer._print(self.k, *args), printer._print(self.args[0], *args), printer._print(self.args[1], *args))

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(2*ImaginaryUnit()*Pi()/2**self.k)]])

    @property
    def name(self):
        return "R%s(%s, %s)" % (self.k, self.args[0], self.args[1])

    @property
    def inputnumber(self):
        return 3

class IRkGate(RkGate):
    """
        Inverse Controlled-Phase Gate
        Does the same thing as the RkGate, but rotates in the opposite direction within the complex plane
        >>> from sympy.physics.qbit import Qbit, IRkGate, apply_gates, QbitZBasisSet
        >>> IRkGate(1,0,2)
        IR2(1, 0)
        >>> from sympy.physics.quantum import represent
        >>> represent(_, QbitZBasisSet(2))
        [1, 0, 0,  0]
        [0, 1, 0,  0]
        [0, 0, 1,  0]
        [0, 0, 0, -I]
        >>> IRkGate(1,0,3)*Qbit(1,1)
        IR3(1, 0)*|11>
        >>> apply_gates(_)
        exp(-pi*I/4)*|11>
    """
    def _apply_ZBasisSet(self, qbits):
        #switch qbit basis and matrix basis when fully implemented
        mat = self.matrix
        args = [self.args[i] for i in reversed(range(2))]
        return self._apply(qbits, mat, args)

    def _sympystr(self, printer, *args):
        return "IR%s(%s, %s)" % (printer._print(self.k, *args), printer._print(self.args[0], *args), printer._print(self.args[1], *args))

    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(-2*ImaginaryUnit()*Pi()/2**self.k)]])

    @property
    def name(self):
        return "IR%s(%s, %s)" % (self.k, self.args[0], self.args[1])

class CTGate(Gate):
    """
        Controlled Pi/8 Gate
        >>> from sympy.physics.qbit import Qbit, CTGate, apply_gates, QbitZBasisSet
        >>> apply_gates(CTGate(0,1)*Qbit(1,1))
        exp(pi*I/4)*|11>
    """
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(I*Pi()/4)]])

class CZGate(Gate):
    """
        Controlled Z-Gate
        >>> from sympy.physics.qbit import Qbit, CZGate, apply_gates, QbitZBasisSet
        >>> apply_gates(CZGate(0,1)*Qbit(1,1))
        -1*|11>
    """
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])


class CPhaseGate(Gate):
    """
        Controlled Phase-Gate
        >>> from sympy.physics.qbit import Qbit, CPhaseGate, apply_gates, QbitZBasisSet
        >>> apply_gates(CPhaseGate(0,1)*Qbit(1,1))
        I*|11>
    """
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,ImaginaryUnit()]])

class CNOTGate(Gate):
    """
        Controlled NOT-Gate (This is the 'entangling gate' most often made refrence to in the literature)
        CNOT, Hadamard and the Pauli-Gates make a universal group in quantum computation
        Flips the second qbit (The target) contingent on the first qbit being 1 (first qbit -> controll bit)

        >>> from sympy.physics.qbit import Qbit, CNOTGate, apply_gates, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(CNOTGate(0,1), QbitZBasisSet(2))
        [1, 0, 0, 0]
        [0, 0, 0, 1]
        [0, 0, 1, 0]
        [0, 1, 0, 0]
        >>> apply_gates(CNOTGate(0,1)*Qbit(0,1))
        |'11'>

    """
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

class SwapGate(Gate):
    """
        SwapGate: swaps two qbits location within the Tensor Product

        >>> from sympy.physics.qbit import Qbit, SwapGate, apply_gates, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(SwapGate(0,1), QbitZBasisSet(2))
        [1, 0, 0, 0]
        [0, 0, 1, 0]
        [0, 1, 0, 0]
        [0, 0, 0, 1]
        >>> apply_gates(SwapGate(0,1)*Qbit(0,1))
        |'10'>
    """
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

class ToffoliGate(Gate):
    """
        ToffoliGate: Also known as the Controlled-Controlled (double controlled) NOT-Gate
        Flips the third qbit (the target) contingent on the first and second qbits being 1 (first && second qbits are controll bits)
        >>> from sympy.physics.qbit import Qbit, ToffoliGate, apply_gates, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(ToffoliGate(1,0,2), QbitZBasisSet(3))
        [1, 0, 0, 0, 0, 0, 0, 0]
        [0, 1, 0, 0, 0, 0, 0, 0]
        [0, 0, 1, 0, 0, 0, 0, 0]
        [0, 0, 0, 0, 0, 0, 0, 1]
        [0, 0, 0, 0, 1, 0, 0, 0]
        [0, 0, 0, 0, 0, 1, 0, 0]
        [0, 0, 0, 0, 0, 0, 1, 0]
        [0, 0, 0, 1, 0, 0, 0, 0]
        >>> apply_gates(ToffoliGate(1,0,2)*Qbit(1,1,1))
        |'011'>
    """
    @property
    def matrix(self):
        return Matrix([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]])

#-----------------------------------------------------------------------------
# Fourier stuff
#-----------------------------------------------------------------------------

class Fourier(Gate):
    """
        Superclass of Quantum Fourier and Inverse Quantum Fourier Gates
        Takes in two args telling which qbits to start and stop doing the (I)QFT on

        >>> from sympy.physics.qbit import Qbit, QFT, IQFT, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(QFT(0,2).decompose(), QbitZBasisSet(2))
        [1/2,  1/2,  1/2,  1/2]
        [1/2,  I/2, -1/2, -I/2]
        [1/2, -1/2,  1/2, -1/2]
        [1/2, -I/2, -1/2,  I/2]
        >>> QFT(0,2).decompose()
        SwapGate(0,1)*HadamardGate(0)*R2(1, 0)*HadamardGate(1)
        >>> represent(IQFT(0,2).decompose(), QbitZBasisSet(2))
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
    def _apply(self, qbits):
        raise NotImplementedError("This command shouldn't happen")

    @property
    def minimumdimension(self):
        #Can apply to a Qbit basisSet up to one less that its last arg
        return self.args[1]-1

    @property
    def inputnumber(self):
        #first input should be start of register, second is end of register
        return 2

    def _represent_ZBasisSet(self, HilbertSize, format = 'sympy'):
        if HilbertSize <= self.minimumdimension:
            raise  HilbertSpaceException("HilbertSize doesn't work")
        product = []
        product.append(eye(2**(self.args[0])))
        product.append(self.matrix)
        product.append(eye(2**(HilbertSize - self.args[1])))
        return TensorProduct(*product)

    @property
    def matrix(self):
        # Can't form the matrix on its own yet
        NotImplementedError("Fourier Transforms don't know how to do this just yet")

    def _apply_ZBasisSet(self, qbits):
        #decomposes self into constituients and applies
        return apply_gates(self.decompose*qbits)

class QFT(Fourier):
    def decompose(self):
        """
            Decomposes QFT into elementary gates
        """
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
        """
            Decomposes IQFT into elementary gates
        """
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
# Miscilaneous Gate sub-classes (Measurement and Controlled Mod)
#-----------------------------------------------------------------------------
class measure(NondistributiveGate):
    """
        A one-shot measurement gate that returns a single outcome from measureing a single gate
    """

    def __new__(cls, *args):
        if len(args) != 1:
            raise QuantumError("Can only measure one at a time for now")
        return Expr.__new__(cls, *args, **{'commutative': False})

    def measure(self, state):
        import random
        qbit = self.args[0]
        state = state.expand()
        prob1 = 0
        #Go through each item in the add and grab its probability
        #This will be used to determine probability of getting a 1
        if isinstance(state, (Mul, Qbit)):
            return state
        for item in state.args:
            if item.args[-1][qbit] == 1:
                prob1 += Mul(*item.args[0:-1])*Mul(*item.args[0:-1]).conjugate()
        if prob1 < random.random():
            choice = 0
        else:
            choice = 1
        result = 0
        for item in state.args:
            if item.args[-1][qbit] == choice:
                result = result + item
        if choice:
            result = result/sqrt(prob1)
        else:
            result = result/sqrt(1-prob1)
        #TODO FIXME I evalf'd becuase of wierdness, need to get measure working in apply for real
        return result.expand().evalf()

class SetZero(NondistributiveGate):
    """
        Gate measures the state of a Qbit within a set of qbits:
            - If the measurement returns 0, it does nothing
            - If the measurement return 1, it applies an XGate
    """
    def measure(self, circuit):
        # Turns item self.args[0]th bit to a zero
        item = self.args[0]
        circuit = measure(item)*circuit
        circuit = apply_gates(circuit)
        circuit = circuit.expand()
        if isinstance(circuit, QAdd):
            if circuit.args[-1][-1][item] == 0:
                return circuit
            part = 0
            for i in circuit.args:
                part = part + QMul(*i.args[:-1])*i.args[-1].flip(item)
            return part
        if isinstance(circuit, Qbit):
            if circuit[item] == 0:
                return circuit
            return circuit.flip(item)
        if isinstance(circuit, QMul):
            if circuit.args[-1][item] == 0:
                return circuit
            return QMul(*circuit.args[:-1])*circuit.args[-1].flip(item)
        return circuit
#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def representHilbertSpace(gateMatrix, HilbertSize, qbits, format='sympy'):
    """
        This is a helper function used by Gates represent functions to represent their matricies in a given HilbertSpace
        e.g. a Hadamard Gate matrix looks different if applied to a different number of gates
    """
    def _operator(first, second, format):
        """
            Returns the Outer product of a one or zero ket and bra
            This is a helper function
        """
        if (first != 1 and first != 0) or (second != 1 and second != 0):
            raise QuantumError("can only make matricies |0><0|, |1><1|, |0><1|, or |1><0|")
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

    def npTensorProduct(*product):
        """
            Wrapper Function that abstracts away numpy kron function as a TensorProduct function
        """
        answer = product[0]
        for item in product[1:]:
            answer = np.kron(answer, item)
        return answer

    if format == 'sympy':
        eye = getattr(gateMatrix, 'eye')
        kron = TensorProduct
    elif format=='numpy':
        #if user specified numpy as matrix format, try to import
        try:
            import numpy as np
        except Exception:
            #If we couldn't load, just revert to sympy
            representHilbertSpace(gateMatrix, HilbertSize, qbits, format ='sympy')
        #redirect eye to the numpy eye function, and kron to a modified numpy function
        gateMatrix = np.matrix(gateMatrix.tolist())
        aeye = getattr(np, 'eye')
        eye = lambda x: np.matrix(aeye(x))
        kron = npTensorProduct
    else:
        raise ValueError()

    if gateMatrix.shape[1] == 2:
        product = []
        qbit = qbits[0]
        #fill product with [I1,Gate,I2] such that the unitaries, I, cause the gate to be applied to the correct qbit
        if qbit != HilbertSize-1:
            product.append(eye(2**(HilbertSize-qbit-1)))
        product.append(gateMatrix)
        if qbit != 0:
            product.append(eye(2**qbit))

        #do the tensor product of these I's and gates
        if format == 'sympy' or format == 'numpy':
            MatrixRep = kron(*product)
        else:
            raise ValueError()
        return MatrixRep

    #If we are dealing with a matrix that is inheritely multi-qubit, do more work
    else:
        #find the control and target qbit(s)
        controls = qbits[:-1]
        controls = [x for x in reversed(controls)]
        target =  qbits[-1]
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
            for j in range(HilbertSize):
                product.append(eye(2))
            n = 0
            #put Operators |0><0|, |1><1|, |0><1|, or |1><0| in place of I's for all control bits
            for item in controls:
                product.pop(HilbertSize-1-item)
                #Operator is picked so that if i = 0xyab (base 2; x,y,a,b = 0 or 1),
                #then operator is |xy><ab| = |x><a|X|y><b| each of (|x><a|) which goes in control bit location
                product.insert(HilbertSize-1-item, _operator(i>>(n+len(controls))&1,(i>>n)&1, format))
                n = n+1
            #put the correct submatrix from matrixarray into target-bit location
            product.pop(HilbertSize-1-target)
            product.insert(HilbertSize-1-target, matrixArray[i])

            #preform Tensor product first time
            if isinstance(answer, (int, Integer)):
                if format == 'sympy' or format == 'numpy':
                    answer = kron(*product)
                else:
                    raise ValueError()
            #add last answer to TensorProduct of what we have
            else:
                if format == 'sympy' or format == 'numpy':
                    answer = answer + kron(*product)
                else:
                    raise ValueError()
        return answer

def TensorProduct(*args):
    """
        Does the Kronecker product of the matricies in the *args tuple (e.g. the tensor product of the matricies)
        For info on how to do Kronecker products, see:
            http://en.wikipedia.org/wiki/Kronecker_product

        >>> from sympy.physics.qbit import TensorProduct
        >>> from sympy.matrices.matrices import Matrix
        >>> TensorProduct(Matrix([0,1]), Matrix([[1,0]]))
        [0, 0]
        [1, 0]
        >>> TensorProduct(Matrix([0,1]), Matrix([[1,0], [0,1]]))
        [0, 0]
        [0, 0]
        [1, 0]
        [0, 1]
    """
    #pull out the first element in the product
    MatrixExpansion  = args[len(args)-1]

    #do the tensor product working from right to left
    for gate in reversed(args[:len(args)-1]):
        rows = gate.rows
        cols = gate.cols
        #go through each row appending tensor product to running MatrixExpansion
        for i in range(rows):
            Start = MatrixExpansion*gate[i*cols]
            #go through each column joining each item
            for j in range(cols-1):
                Start = Start.row_join(MatrixExpansion*gate[i*cols+j+1])
            #if this is the first element in row, make it the start of the new row
            if i == 0:
                Next = Start
            else:
                Next = Next.col_join(Start)
        MatrixExpansion = Next
    return MatrixExpansion

def apply_gates(circuit, basis = QbitZBasisSet(1), floatingPoint = False):
    """
        Uses the Gate sequence to map initial state into final state without creating a representative matrix first
        Thus, It maps HadamardGate(0)*Qbit(0,0) -> Qbit(0,1)/sqrt(2) + Qbit(0,0)/sqrt(2)

        >>> from sympy.physics.qbit import HadamardGate, Qbit, apply_gates
        >>> apply_gates(HadamardGate(0)*Qbit(0,1))
        -2**(1/2)/2*|01> + 2**(1/2)/2*|00>
    """
    # if all we have is a Qbit without any gates, return the qbit
    if isinstance(circuit, Qbit):
        return circuit

    #if we have a Mul object, get the state of the system
    elif isinstance(circuit, QMul):
        states = circuit.args[len(circuit.args)-1]
        states = states.expand()

    #if we have an add object with gates mixed in, apply_gates recursively
    elif isinstance(circuit, QAdd):
        result = 0
        for i in circuit.args:
            result = result + apply_gates(i, basis, floatingPoint)
        if floatingPoint:
            result = result.evalf()
        return result

    state_coeff = 1
    #pick out each object that multiplies the state
    for multiplier in reversed(circuit.args[:len(circuit.args)-1]):

        #if the object that mutliplies is a Gate, we will apply it once
        if isinstance(multiplier, Gate):
            gate = multiplier
            number_of_applications = 1

        #if the object that multiplies is a Pow who's base is a Gate, we will apply Pow.exp times
        elif isinstance(multiplier, QPow) and isinstance(multiplier.base, Gate):
            gate = multiplier.base
            number_of_applications = multiplier.exp

        #if the object that multiplies is not a gate of any sort, we apply it by multiplying
        else:
            state_coeff = multiplier*state_coeff
            continue

        #if states is in superposition of states (a sum of qbits states), applyGates to each state contined within
        if isinstance(states, QAdd):
            #special check for non-distributivity, do all at once
            if isinstance(gate, NondistributiveGate):
                states = gate.measure(states)
            else:
                result = 0
                for state in states.args:
                    result = result + apply_gates(gate**number_of_applications*state, basis, floatingPoint)
                if floatingPoint:
                    result = result.evalf()
                states = result
                states = states.expand()

        #if we have a mul, apply gate to each register and multiply result
        elif isinstance(states, QMul):
            #find the Qbits in the Mul
            for i in range(len(states.args)):
                if isinstance(states.args[i],Qbit):
                    break
            #if we didn't find one, something is wrong
            if not isinstance(states.args[i],Qbit):
                print states
                raise QuantumError()

            #apply the gate the right number of times to this state
            coefficient = Mul(*(states.args[:i]+states.args[i+1:]))#TODO
            states = apply_gates(gate**(number_of_applications)*states.args[i], basis, floatingPoint)
            states = coefficient*states
            states = states.expand()

        #If we have a single Qbit, apply to this Qbit
        elif isinstance(states, Qbit):
            if isinstance(gate, NondistributiveGate):
                states = gate.measure(states)
            else:
                basis_name = basis.__class__.__name__
                apply_method_name = '_apply_%s' % basis_name
                apply_method = getattr(gate, apply_method_name)
                states = apply_method(states)
                states = states.expand()
                number_of_applications -= 1
                while number_of_applications > 0:
                    states = apply_gates(gate*states, basis, floatingPoint)
                    number_of_applications -= 1

        #if it's not one of those, there is something wrong
        else:
            raise QuantumError()

    #tack on any coefficients that were there before and simplify
    states = state_coeff*states
    if isinstance(states, (Add,Pow, Mul)):
        states = states.expand()
    return states

def matrix_to_qbits(matrix):
    """
        Takes a matrix representation of a qbit and puts in into a sum of Qbit eigenstates of the ZBasisSet
        >>> from sympy.physics.qbit import matrix_to_qbits, Qbit, QbitZBasisSet
        >>> from sympy.physics.quantum import represent
        >>> represent(Qbit(0,1), QbitZBasisSet(2))
        [0]
        [1]
        [0]
        [0]
        >>> matrix_to_qbits(represent(Qbit(0,1), QbitZBasisSet(2)))
        |'01'>
    """
    #make sure it is of correct dimensions for a qbit-matrix representation
    qbit_number = log(matrix.rows,2)
    if matrix.cols != 1 or not isinstance(qbit_number, Integer):
        raise QuantumError()

    #go through each item in matrix, if element is not zero, make it into a qbit item times coefficient
    result = 0
    mlistlen = len(matrix.tolist())
    for i in range(mlistlen):
        if matrix[i] != 0:
            #form qbit array; 0 in bit-locations where i is 0, 1 in bit-locations where i is 1
            qbit_array = [1 if i&(1<<x) else 0 for x in range(qbit_number)]
            qbit_array.reverse()
            result = result + matrix[i]*Qbit(*qbit_array)

    #if sympy simplified by pulling out a constant coefficient, undo that
    if isinstance(result, (Mul,Add,Pow)):
        result = result.expand()
    return result

def qbits_to_matrix(qbits):
    """
        Coverts a QAdd of qbit objects into

        >>> from sympy.physics.qbit import qbits_to_matrix, Qbit
        >>> from sympy.functions.elementary.miscellaneous import sqrt
        >>> qbits_to_matrix(Qbit(0,0)/sqrt(2) + Qbit(0,1)/sqrt(2))
        [2**(1/2)/2]
        [2**(1/2)/2]
        [         0]
        [         0]
    """
    #get rid of multiplicative constants
    qbits = qbits.expand()

    #if we have a Mul object, find the qbit part qbits to matrix it
    if isinstance(qbits, QMul):
        for i in range(len(qbits.args)):
            if isinstance(qbits.args[i], Qbit):
                break
        if not isinstance(qbits.args[i], Qbit):
            raise QuantumError()
        #recursively turn qbit into matrix
        return Mul(*(qbits.args[:i] + qbits.args[i+1:]))*qbits_to_matrix(qbits.args[i])
    #recursively turn each item in an add into a matrix
    elif isinstance(qbits, QAdd):
        result = qbits_to_matrix(qbits.args[0])
        for element in qbits.args[1:]:
            result = result + qbits_to_matrix(element)
        return result
    #if we are at the bottom of the recursion, have the base case be representing the matrix
    elif isinstance(qbits, Qbit):
        return qbits._represent_QbitZBasisSet(QbitZBasisSet(len(qbits))) #TODO other bases with getattr
    else:
        raise QuantumError("Malformed input")

def gatesimp(circuit):
    """
        Simplifies gates symbolically
        >>> from sympy.physics.qbit import gatesimp, HadamardGate
        >>> gatesimp(HadamardGate(1)**3*HadamardGate(0))
        HadamardGate(0)*HadamardGate(1)
    """

    #bubble sort(?) out gates that commute
    circuit = gatesort(circuit)

    #do simplifications
    if isinstance(circuit, QMul):
        for i in range(len(circuit.args)):
            #H,X,Y or Z squared is 1. T**2 = S, S**2 = Z
            if isinstance(circuit.args[i], QPow):
                if isinstance(circuit.args[i].base, (HadamardGate, XGate, YGate, ZGate)) and isinstance(circuit.args[i].exp, Number):
                    newargs = (circuit.args[:i] + (circuit.args[i].base**(circuit.args[i].exp % 2),) + circuit.args[i+1:])
                    circuit = gatesimp(QMul(*newargs))
                    break
                elif isinstance(circuit.args[i].base, PhaseGate):
                    newargs = (circuit.args[:i] + (ZGate(circuit.args[i].base.args[0])**(Integer(circuit.args[i].exp/2)), circuit.args[i].base**(circuit.args[i].exp % 2)) + circuit.args[i+1:])
                    circuit =  gatesimp(QMul(*newargs))
                    break
                elif isinstance(circuit.args[i].base, TGate):
                    newargs = (circuit.args[:i] + (SGate(circuit.args[i].base.args[0])**Integer(circuit.args[i].exp/2), circuit.args[i].base**(circuit.args[i].exp % 2)) + circuit.args[i+1:])
                    circuit =  gatesimp(QMul(*newargs))
                    break

    return circuit

def gatesort(circuit):
    """
        Sorts the gates while keeping track of commutation relations
        (things that apply to the same qbit do not commute with each other)

        >>> from sympy.physics.qbit import HadamardGate, XGate, YGate, CNOTGate, gatesort
        >>> gatesort(YGate(2)**2*HadamardGate(0)*CNOTGate(0,1)*XGate(1)*YGate(0))
        HadamardGate(0)*CNOTGate(0,1)*YGate(0)*XGate(1)*(YGate(2))**2
    """
    #bubble sort of gates checking for commutivity of neighbor (Python doesn't have a do-while)
    changes = True
    while changes:
        changes = False
        cirArray = circuit.args
        for i in range(len(cirArray)-1):
            #Go through each element and switch ones that are in wrong order
            if isinstance(cirArray[i], (Gate, QPow)) and isinstance(cirArray[i+1], (Gate, QPow)):
                if isinstance(cirArray[i], QPow):
                    first = cirArray[i].base
                else:
                    first = cirArray[i]

                if isinstance(cirArray[i+1], QPow):
                    second = cirArray[i+1].base
                else:
                    second = cirArray[i+1]

                if first.args > second.args:
                    #make sure elements commute, meaning they do not affect ANY of the same qbits
                    commute = True
                    for arg1 in cirArray[i].args:
                       for arg2 in cirArray[i+1].args:
                            if arg1 == arg2:
                                commute = False
                    # if they do commute, switch them
                    if commute:
                        circuit = QMul(*(circuit.args[:i] + (circuit.args[i+1],) + (circuit.args[i],) + circuit.args[i+2:]))
                        cirArray = circuit.args
                        changes = True
    return circuit

#-----------------------------------------------------------------------------
# Quantum ALU contains quantum versions of arithmetic operations
# These use elementary gates to preform arithmetic operations
#-----------------------------------------------------------------------------

def ADD(InReg, InOutReg, carryReg):
    """
        Add's the InReg to the InOutReg, and stores the result in InOutReg
        CarryReg is used as extra storage space to preform calculation

        Doctest does 2 + 2 == 4:
        >>> from sympy.physics.qbit import ADD, apply_gates, ToffoliGate, SetZero, CNOTGate, Qbit
        >>> apply_gates(ADD((0,1,2,3), (4,5,6,7), (8,9,10,11))*Qbit(1,1,1,1,0,0,1,0,0,0,1,0))
        |'000001000010'>
    """
    circuit = 1
    for item in carryReg:
        circuit = SetZero(item)*circuit
    for i in range(len(InReg)):
        if i != (len(InReg)-1):
            circuit = ToffoliGate(InReg[i], InOutReg[i], carryReg[i+1])*circuit
            circuit = ToffoliGate(InReg[i], carryReg[i], carryReg[i+1])*circuit
            circuit = ToffoliGate(InOutReg[i], carryReg[i], carryReg[i+1])*circuit
        circuit = CNOTGate(InReg[i], InOutReg[i])*circuit
        circuit = CNOTGate(carryReg[i], InOutReg[i])*circuit
    for item in carryReg:
        circuit = SetZero(item)*circuit
    return circuit

def ControlledADD(InReg, InOutReg, carryReg, control):
    """
        Same as ADD, but only does the add if control bit is true
        Add's InReg and InOutReg to each other and stores result in InOutReg contingent on control being true
        >>> from sympy.physics.qbit import ControlledADD, Qbit, apply_gates
        >>> apply_gates(ControlledADD((0,1,2,3), (4,5,6,7), (8,9,10,11), 12)*Qbit(1,0,0,0,0,0,0,1,0,0,0,1,0))
        |'1000001000010'>
        >>> apply_gates(ControlledADD((0,1,2,3), (4,5,6,7), (8,9,10,11), 12)*Qbit(0,0,0,0,0,0,0,1,0,0,0,1,0))
        |'0000000100010'>
    """
    circuit = 1
    for item in carryReg:
        circuit = SetZero(item)*circuit
    for i in range(len(InReg)):
        if i != (len(InReg)-1):
            circuit = ToffoliGate(InReg[i], InOutReg[i], carryReg[i+1])*circuit
            circuit = ToffoliGate(InReg[i], carryReg[i], carryReg[i+1])*circuit
            circuit = ToffoliGate(InOutReg[i], carryReg[i], carryReg[i+1])*circuit
        circuit = ToffoliGate(control, InReg[i], InOutReg[i])*circuit
        circuit = ToffoliGate(control, carryReg[i], InOutReg[i])*circuit
    for item in carryReg:
        circuit = SetZero(item)*circuit
    return circuit

def Bitshift(Register, number, tempStorage):
    """
        Shifts the given Register by the given number
        tempStorage is the single qbit of storage needed
        >>> from sympy.physics.qbit import Bitshift, Qbit, apply_gates
        >>> apply_gates(Bitshift((0,1,2,3,4), 1, 5)*Qbit(0,0,0,0,1,0))
        |'000100'>
        >>> apply_gates(Bitshift((0,1,2,3,4), -1, 5)*Qbit(0,0,0,0,1,0))
        |'000001'>
    """
    circuit = SetZero(tempStorage)
    if number > 0:
        for i in range(abs(number)):
            circuit = SwapGate(Register[-1], tempStorage)*circuit
            for i in reversed(range(len(Register)-1)):
                circuit = SwapGate(i, i+1)*circuit
            circuit = SetZero(tempStorage)*circuit
        return circuit
    elif number < 0:
        for i in range(abs(number)):
            circuit = SwapGate(Register[0], tempStorage)*circuit
            for i in range(len(Register)-1):
                circuit = SwapGate(i, i+1)*circuit
            circuit = SetZero(tempStorage)*circuit
        return circuit
    else:
        return 1

def Multiply(InReg1, InReg2, OutReg, carryReg):
    """
        Multiplies the first register by the second and puts the result in OutReg
        CarryReg is used to keep carry information

        Here is 2*2
        >>> from sympy.physics.qbit import Multiply, Qbit, apply_gates
        >>> apply_gates(Multiply((0,1,2,3), (4,5,6,7), (8,9,10,11), (12,13,14,15))*Qbit(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0))
        |'00000010000100000'>
    """
    circuit = 1
    for item in OutReg:
        circuit = SetZero(item)*circuit
    for item in carryReg:
        circuit = SetZero(item)*circuit
    for i in InReg2:
        circuit = ControlledADD(InReg1, OutReg, carryReg, i)*circuit
        circuit = Bitshift(InReg1, 1, carryReg[0])*circuit
    return circuit
