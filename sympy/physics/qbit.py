"""
    Single qbits and their gates
"""
from sympy.physics.hilbert import l2
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
import math


class QbitZBasisSet(BasisSet):
    
    def __init__(self, nqbits):
        self.hilbert_space = l2(2)**nqbits
    
    def __new__(cls, nqbits):
        return BasisSet.__new__(cls, nqbits)

class QbitXBasis(BasisSet):
    
    def __init__(self, nqbits):
        self.hilbert_space = l2(2)**nqbits
    
    def __new__(cls, nqbits):
        return BasisSet.__new__(cls, nqbits)

class QbitYBasis(BasisSet):
    
    def __init__(self, nqbits):
        self.hilbert_space = l2(2)**nqbits
    
    def __new__(cls, nqbits):
        return BasisSet.__new__(cls, nqbits)

class Qbit(Ket):
    """
        Represents a single definite quantum state
    """
    outDecimal = False
    def __new__(cls, *args, **options):
        import math
        #If they just give us one number, express it in the least number of bits possible
        if args[0] > 1 and len(args) == 1:
            array = [(args[0]>>i)&1 for i in reversed(range(int(math.ceil(math.log(args[0], 2)+.01)+.001)))]
            array = sympify(array)
            obj = Expr.__new__(cls, *array)
            return obj
        #if they give us two numbers, the second number is the number of bits on which it is expressed)
        #Thus, Qbit(0,5) == |00000>. second argument can't be one becuase of intersection and uslessesness of possibility
        elif len(args) == 2 and args[1] > 1:
            array = [(args[0]>>i)&1 for i in reversed(range(args[1]))]
            array = sympify(array)
            obj = Expr.__new__(cls, *array)
            return obj
        for element in args:
            if not (element == 1 or element == 0):
                raise Exception("Values must be either one or zero")
        args = sympify(args)
        obj = Expr.__new__(cls, *args)
        return obj

    @property
    def hilbert_space(self):
        return l2(2)
    
    @property
    def dimension(self):
        return len(self.args)
    
    def __len__(self):
        return self.dimension
    
    def __getitem__(self, bit):
        if bit > self.dimension - 1:
            raise Exception()
        return self.args[int(self.dimension-bit-1)]

    def _print_name(self, printer, *args):
        string = self.to_string()
        return printer._print(string, *args)

    def _print_name_pretty(self, printer, *args):
        string = self.to_string()
        pform = printer._print(string, *args)
        return pform

    def to_string(self):
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
        #check now needed TODO
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
        return 1/sqrt(2)*Matrix([[1,1],[1,-1]])
    
    @property
    def YBasisTransform(self):
        return Matrix([[ImaginaryUnit(),0],[0,-ImaginaryUnit()]])

class QbitX(Qbit):
    def __new__(cls, *args):
        for element in args:
            if not (element == '+' or element == '-'):
                raise Exception("Values must be either + or -")
        return Expr.__new__(cls, *args, **{'commutative': False})

class Gate(Operator):
    """
    A gate operator that acts on qubit(s)
    (will need tensor product to get it working for multiple Qbits)
    """
    name = 'Gate'
    def __new__(cls, *args):
        args = sympify(args)
        obj = Expr.__new__(cls, *args)
        if obj.inputnumber != len(args):
            num = obj.inputnumber
            raise Exception("This gate applies to %d qbits" % (num))
        for i in range(len(args)):
            if args[i] in (args[:i] + args[i+1:]):
                raise Exception("Can't have duplicate control and target bits!")
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
            raise Exception("can't apply to object that is not a qbit")
        return result
    
    def _apply_QbitZBasisSet(self, qbits):
        #switch qbit basis and matrix basis when fully implemented
        mat = self.matrix
        args = [self.args[i] for i in reversed(range(len(self.args)))]
        return self._apply(qbits, mat, args)
    
    def _sympyrepr(self, printer, *args):
        return "%s(%s)" %  (printer._print(self.__class__.__name__, *args), printer._print(self.args, *args))
    
    def _represent_QbitZBasisSet(self, basis, format = 'sympy'):
        if self.minimumdimension >= basis.dimension:
            raise HilbertSpaceException()
        gate = self.matrix
        if basis.dimension  == 1:
            return gate
        else:
            m = representHilbertSpace(gate, basis.dimension, self.args, format)
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
        return 1/sqrt(2)*Matrix([[1,1],[1,-1]])
    
    @property
    def YBasisTransform(self):
        return Matrix([[ImaginaryUnit(),0],[0,-ImaginaryUnit()]])

#class keeps track of fact that it can't be distributed
class NondistributiveGate(Gate):
    def __new__(cls, *args):
        args = sympify(args)
        if len(args) != 1:
            raise Exception("Can only measure one at a time for now")
        return Expr.__new__(cls, *args, **{'commutative': False})
    
    def measure(states):
        raise NotImplementedError("This type of measure has not been implemented")

def representHilbertSpace(gateMatrix, HilbertSize, qbits, format='sympy'):
    
    #returns |first><second|, if first and second are 0 or 1.
    def _operator(first, second, format):
        if (first != 1 and first != 0) or (second != 1 and second != 0):
            raise Exception("can only make matricies |0><0|, |1><1|, |0><1|, or |1><0|")
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
    
    #Wrapper function that gives np.kron same interface as TensorProduct
    def npTensorProduct(*product):
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

#This is black box, I need to find a way to do this more realistically FIXME TODO
class controlledMod(Gate):
    def __new__(cls, *args):
        return Expr.__new__(cls, *args, **{'commutative': False})
    
    def _apply_QbitZBasisSet(self, qbits):
        t = self.args[0]
        a = self.args[1]
        N = self.args[2]
        n = 1
        k = 0
        for i in range(t):
            k = k + n*qbits[t+i]
            n = n*2
        out = int(a**k%N)
        outarray = list(qbits.args[0:t])
        for i in reversed(range(t)):            
            outarray.append((out>>i)&1)
        return Qbit(*outarray)

class ensembleMeasure(NondistributiveGate): #for now, ensemble measure does not do partial measurements FIXIT
    def __new__(cls, *args):
        return Expr.__new__(cls, *args, **{'commutative': False})
    
    def measure(self, state):
        state = state.expand()
        if isinstance(state, (Qbit, Mul)):
                return state
        if len(self.args) == 0:
            #Do whole measuremnet
            retVal = []
            for each_item in state.args:
                retVal.append((Mul(*each_item.args[0:-1])*Mul(*each_item.args[0:-1]).conjugate(), each_item.args[-1]))
            return retVal
        else:
            raise NotImplementedError("Don't have partial ensemble measures done yet, check back later")
            #Do partial measurement

class measure(NondistributiveGate):
    def __new__(cls, *args):
        if len(args) != 1:
            raise Exception("Can only measure one at a time for now")
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

class RkGate(Gate):
    def __new__(cls, *args):
        obj = Expr.__new__(cls, *args, **{'commutative': False})
        if obj.inputnumber+1 != len(args):
            num = obj.inputnumber
            raise Exception("This gate applies to %d qbits" % (num))
        return obj
    
    def _apply_QbitZBasisSet(self, qbits):
        #switch qbit basis and matrix basis when fully implemented
        mat = Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(2*ImaginaryUnit()*Pi()/2**self.args[2])]])
        args = [self.args[i] for i in reversed(range(2))]
        return self._apply(qbits, mat, args)
    
    def _sympystr(self, printer, *args):
        return "R%s(%s, %s)" % (printer._print(self.args[2], *args), printer._print(self.args[0], *args), printer._print(self.args[1], *args))
    
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(2*ImaginaryUnit()*Pi()/4)]])
    
    def _represent_ZBasisSet(self, HilbertSize, format = 'sympy'):
        if self.minimumdimension >= HilbertSize:
            raise HilbertSpaceException()
        gate = Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(2*ImaginaryUnit()*Pi()/2**self.args[2])]])
        if HilbertSize  == 1:
            return gate
        else:
            m = representHilbertSpace(gate, HilbertSize, self.args[:2], format)
            return m
    
    @property
    def minimumdimension(self):
        return max(self.args[:2])

class IRkGate(RkGate):
    def _apply_ZBasisSet(self, qbits):
        #switch qbit basis and matrix basis when fully implemented
        mat = Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(-2*ImaginaryUnit()*Pi()/2**self.args[2])]])
        args = [self.args[i] for i in reversed(range(2))]
        return self._apply(qbits, mat, args)
    
    def _sympystr(self, printer, *args):
        return "IR%s(%s, %s)" % (printer._print(self.args[2], *args), printer._print(self.args[0], *args), printer._print(self.args[1], *args))
    
    def _represent_ZBasisSet(self, HilbertSize, format = 'sympy'):
        if self.minimumdimension >= HilbertSize:
            raise HilbertSpaceException()
        gate = Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(-2*ImaginaryUnit()*Pi()/2**self.args[2])]])
        if HilbertSize  == 1:
            return gate
        else:
            m = representHilbertSpace(gate, HilbertSize, self.args[:2], format)
            return m

class CZGate(Gate):
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])

class SwapGate(Gate):
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

class CPhaseGate(Gate):
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,ImaginaryUnit()]])

class ToffoliGate(Gate):
    @property
    def matrix(self):
        return Matrix([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]])

class CNOTGate(Gate):
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

class HadamardGate(Gate):
    """
    An object representing a Hadamard Gate
    """
  
    @property
    def matrix(self):
        from sympy.functions.elementary.miscellaneous import sqrt
        return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))

class XGate(Gate):
    """
    An object representing a Pauli-X gate:
    """
    @property
    def matrix(self):
        return Matrix([[0, 1], [1, 0]])

class YGate(Gate):
    """
    An object representing a Pauli-Y gate:
    """
   
    @property
    def matrix(self):
        return Matrix([[0, complex(0,-1)], [complex(0,1), 0]])

class CTGate(Gate):
    @property
    def matrix(self):
        return Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,exp(I*Pi()/4)]])

class ZGate(Gate):
    """
    An object representing a Pauli-Z gate:
    """
   
    @property
    def matrix(self):
        return Matrix([[1, 0], [0, -1]])

class PhaseGate(Gate):
    """
    An object representing a phase gate:
    """
   
    @property
    def matrix(self):
        return Matrix([[1, 0], [0, complex(0,1)]])

class TGate(Gate):
    """
    An object representing a pi/8 gate:
    """
    
    @property
    def matrix(self):
        return Matrix([[1, 0], [0, exp(I*Pi()/4)]])

def apply_gates(circuit, basis = QbitZBasisSet(1), floatingPoint = False):
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
                raise Exception()
            
            #apply the gate the right number of times to this state
            coefficient = states._new_rawargs(states.evaluates, states.hilbert_space, *(states.args[:i]+states.args[i+1:]))#TODO
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
            raise Exception()
    
    #tack on any coefficients that were there before and simplify
    states = state_coeff*states
    if isinstance(states, (Add,Pow, Mul)):
        states = states.expand()
    return states
    

def matrix_to_qbits(matrix):
    #make sure it is of correct dimensions for a qbit-matrix representation
    qbit_number = log(matrix.rows,2)
    if matrix.cols != 1 or not isinstance(qbit_number, Integer):
        raise Exception()
    
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
    #get rid of multiplicative constants
    qbits = qbits.expand()
    
    #if we have a Mul object, find the qbit part qbits to matrix it
    if isinstance(qbits, QMul):
        for i in range(len(qbits.args)):
            if isinstance(qbits.args[i], Qbit):
                break
        if not isinstance(qbits.args[i], Qbit):
            raise Exception()
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
        raise Exception("Malformed input")

def gatesimp(circuit):
    """ will simplify gates symbolically"""
    #Pull gates out of inner Add's and Mul's?
    
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
            #Deal with HXH = Z, HZH = X, HYH = -Y
            if isinstance(circuit.args[i], HadamardGate):
                #check for X,Y,Z in front
                pass
                #check for Hadamard to right of that
                #replace stuff
    
    return circuit

def gatesort(circuit):
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

class HilbertSpaceException(Exception):
    pass

class Fourier(Gate):
    def __new__(self, *args):
        if args[0] >= args[1]:
            raise Exception("Start must be smaller than finish")
        return Expr.__new__(self, *args)
    
    @property
    def _apply(self, qbits):
        raise NotImplementedError("This command doesn't make sense for a Fourier Transform")
    
    @property
    def minimumdimension(self):
        return self.args[1]-1
    
    @property
    def inputnumber(self):
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
        N = 2**(self.args[1]-self.args[0])
        if isinstance(self, QFT):
            omega = exp(2*I*Pi()/(N))
        else:
            omega =  exp(-2*I*Pi()/(N))
        mat = 0
        for i in range(N):
            temp = 0
            for j in range(N):
                if temp == 0:
                    temp = Matrix([[1]])
                else:
                    temp = temp.row_join(Matrix([[omega**(i*j%N)]]))
            if mat == 0:
                mat = temp
            else:
                mat = mat.col_join(temp)
        mat = mat/sqrt(N)
        return mat
    
    def _apply_ZBasisSet(self, qbits):
        mat = qbits_to_matrix(qbits)
        return self.DFT(mat)
    
    def DFT(self, qbitmat):
        N = qbitmat.rows
        if isinstance(self, QFT):
            omega = exp(2*I*Pi()/(N))
        else:
            omega = exp(-2*I*Pi()/(N))
        retVal = 0
        for i in range(N):
            if qbitmat[i] == 0:
                continue
            temp = [omega**(j*i)/sqrt(N) for j in range(N)]
            if retVal == 0:
                retVal = Matrix(temp)
            else:
                retVal = retVal + Matrix(temp)
        return retVal

####### ALU #########

class ADD(Gate): #TODO check what happens when we carry for diff number of in-out sizes
    def __new__(cls, *args):
        InReg = args[0]
        InOutReg = args[1]
        carryReg = args[2]
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

class ControlledADD(Gate): #TODO check what happens when we carry for diff number of in-out sizes
    def __new__(cls, *args):
        InReg = args[0]
        InOutReg = args[1]
        carryReg = args[2]
        control = args[3]
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

class Bitshift(Gate): #check
    def __new__(cls, *args):
        Register = args[0]
        number = args[1]
        tempStorage = args[2]
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

class Multiply(NondistributiveGate): #This is harder than I thought will need to fix-it (Controlled-Add?)
    def __new__(cls, *args):
        if len(args) != 4:
            raise Exception("Must input InReg1, InReg2, OutReg, carryReg Tuples")
        InReg1 = args[0]
        InReg2 = args[1]
        OutReg = args[2]
        carryReg = args[3]
        circuit = 1
        for item in OutReg:
            circuit = SetZero(item)*circuit
        for item in carryReg:
            circuit = SetZero(item)*circuit
        for i in InReg2:
            circuit = ControlledADD(InReg1, OutReg, carryReg, i)*circuit
            circuit = Bitshift(InReg1, 1, carryReg[0])*circuit
        return circuit

class SetZero(NondistributiveGate):
    def measure(self, circuit):
        item = self.args[0]
        circuit = measure(item)*circuit
        circuit = apply_gates(circuit)
        circuit = circuit.expand()
        if isinstance(circuit, Add):
            if circuit.args[-1].args[-1][item] == 0:
                return circuit
            part = 0
            for i in circuit.args:
                part = part + i.args[:-1]*i.args[-1].flip(item)
            return part
        if isinstance(circuit, Qbit):
            if circuit[item] == 0:
                return circuit
            return circuit.flip(item)
        if isinstance(circuit, Mul):
            if circuit.args[-1][item] == 0:
                return circuit
            return circuit.args[:-1]*circuit.args[-1].flip(item)
        return circuit

class QFT(Fourier):
    def decompose(self):
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
         
