"""
    Single qbits and their gates
"""
from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify, Matrix, elementary
from sympy.core.numbers import *
from sympy.core.basic import S, sympify
from sympy.core.function import Function
from sympy.functions.elementary.exponential import *
from sympy.functions.elementary.miscellaneous import *
from sympy.matrices.matrices import *
from sympy.simplify import *
import copy

class Qbit(Expr):
    """
    Represents a single quantum gate    
    """
    def __new__(cls, *args):
        obj = Expr.__new__(cls, *args, commutative = False)
        return obj

    @property
    def dimension(self):
        return len(self.args)

    def __len__(self):
        return self.dimension

    def __getitem__(self, bit):
        if bit > self.dimension - 1:
            raise Exception()
        return self.args[self.dimension-bit-1]

    def _sympystr(self, printer, *args):
        string = ""
        for it in self.args:
            string = string + str(it)
        return "|%s>" % printer._print(string, *args)

    def _sympyrepr(self, printer, *args):
        return "%s(%i)" %  (self.__class__.__name__, printer._print(self.name, args))

    def flip(self, *args):
        newargs = list(self.args[:])
        for i in args:
            if newargs[self.dimension-i-1] == 1:
                newargs[self.dimension-i-1] = 0
            else:
                newargs[self.dimension-i-1] = 1
        return Qbit(*newargs)

    def _represent(self):
        n = 1
        definiteState = 0
        for it in reversed(self.args):
            definiteState += n*it
            n = n*2
        result = [0 for x in range(2**self.dimension)]
        result[definiteState] = 1

        return Matrix(result) 

class Gate(Expr):
    """
    A gate operator that acts on qubit(s)
    (will need tensor product to get it working for multiple Qbits)
    """
    def __new__(cls, *args):
        obj = Expr.__new__(cls, *args, commutative = False)
        return obj

    @property
    def matrix(self):
        raise NotImplementedError("matrixRep Not implemented")
        
    @property
    def minimumdimension(self):
        return max(self.args)

    def _apply_ZBasisSet(self, qbits):
        assert isinstance(qbits, Qbit), "msg"
        # check number of qbits this gate acts on
        # check to make sure dimensions are OK
        target_bit = self.args[0]
        mat = self.matrix
        if isinstance(qbits, Qbit):
            column_index = qbits[target_bit]
            column = mat[:,column_index]
            if column_index:
                qbits = column[1]*qbits + column[0]*qbits.flip(target_bit)
            else:
                qbits = column[0]*qbits + column[1]*qbits.flip(target_bit)
        else:
            raise Exception()
        return qbits 

    def _sympyrepr(self, printer, *args):
        return "%s(%s)" %  (self.__class__.__name__, printer._print(self.name, args))

    def _represent_ZBasisSet(self, HilbertSize):
        if self.minimumdimension >= HilbertSize:
            raise HilbertSpaceException()
        gate = self.matrix
        if HilbertSize  == 1:            
            return gate
        else:
            m = representHilbertSpace(gate, HilbertSize, self.args[0])
            return m 

    def _represent_XBasisSet(self, HilbertSize):
        raise NotImplementedError("X-Basis Representation not implemented")

    def _represent_YBasisSet(self, HilbertSize):
        raise NotImplementedError("Y-Basis Representation not implemented")

def representHilbertSpace(gateMatrix, HilbertSize, qbit, format='sympy'):
    """   if format=='sympy':
    elif format=='numpy':
    else:
        raise ValueError()
    """
    product = []
    #fill product with [I1,Gate,I2] such that the unitaries, I, cause the gate to be applied to the correct qbit  
    if qbit != HilbertSize-1:
        product.append(eye(2**(HilbertSize-qbit-1)))    
    product.append(gateMatrix)
    if qbit != 0:
        product.append(eye(2**qbit))

    #do the tensor product of these I's and gates
    MatrixRep = TensorProduct(*product)
    return MatrixRep

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

class BasisSet(Expr):
    pass

class ZBasisSet(BasisSet):
    pass

class XBasisSet(BasisSet):
    pass

class YBasisSet(BasisSet):
    pass

class HadamardGate(Gate):
    """
    An object representing a Hadamard Gate
    """
    def _sympystr(self, printer, *args):
        return "H(%s)" % printer._print(self.args[0], *args)

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))

    def _represent_YBasisSet(self, HilbertSize):
        pass

    @property
    def matrix(self):
        return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))
                  
class XGate(Gate):
    """
    An object representing a Pauli-X gate:
    """
    def _sympystr(self, printer, *args):
        return "X(%s)" % printer._print(self.args[0], *args)

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[1,0],[0,-1]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

    @property
    def matrix(self):
        return Matrix([[0, 1], [1, 0]])

class YGate(Gate):
    """
    An object representing a Pauli-Y gate:
    """
    def _sympystr(self, printer, *args):
        return "Y(%s)" % printer._print(self.args[0], *args)

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[0,complex(0,1)],[complex(0,-1),0]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

    @property
    def matrix(self):
        return Matrix([[0, complex(0,-1)], [complex(0,1), 0]])

class ZGate(Gate):
    """
    An object representing a Pauli-Z gate:
    """
    def _sympystr(self, printer, *args):
        return "Z(%s)" % printer._print(self.args[0], *args)

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[0,1],[1,0]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, -1]])

class PhaseGate(Gate):
    """
    An object representing a phase gate:
    """
    def _sympystr(self, printer, *args):
        return "S(%s)" % printer._print(self.args[0], *args)

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[complex(.5,.5), complex(.5,-.5)], [complex(.5,-.5),complex(.5,.5)]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

    @property
    def matrix(self):
        return Matrix([[1, 0], [0, complex(0,1)]])

class TGate(Gate):
    """
    An object representing a pi/8 gate:
    """
    def _sympystr(self, printer, *args):
        return "T(%s)" % printer._print(self.args[0], *args)

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[.5+.5*exp(complex(0,Pi/4)),.5-.5*exp(complex(0,Pi/4))],[.5-.5*exp(complex(0,Pi/4)),.5+.5*exp(complex(0,Pi/4))]])

    def _represent_YBasisSet(self, HilbertSize):
        pass


    @property
    def matrix(self):
        return Matrix([[1, 0], [0, exp(I*Pi()/4)]])

def apply_gates(circuit, basis = ZBasisSet()):
    # if all we have is a Qbit without any gates, return
    if isinstance(circuit, Qbit):
        return circuit

    #if we have a Mul object, get the state of the system
    if isinstance(circuit, Mul):
        states = circuit.args[len(circuit.args)-1]
        states = states.expand()

    #if we have an add object with gates mixed in, apply_gates recursively
    if isinstance(circuit, Add):
        result = 0
        for i in circuit.args:
            result = result + apply_gates(i, basis)
        return result        

    state_coeff = 1
    #pick out each object that multiplies the state
    for multiplier in reversed(circuit.args[:len(circuit.args)-1]):
    
        #if the object that mutliplies is a Gate, we will apply it once
        if isinstance(multiplier, Gate):
            gate = multiplier
            number_of_applications = 1
            
        #if the object that multiplies is a Pow who's base is a Gate, we will apply Pow.exp times
        elif isinstance(multiplier, Pow) and isinstance(multiplier.base, Gate):
            gate = multiplier.base
            number_of_applications = multiplier.exp
            
        #if the object that multiplies is not a gate of any sort, we apply it by multiplying
        else:
            state_coeff = multiplier*state_coeff
            continue

        #if states is in superposition of states (a sum of qbits states), applyGates to each state contined within
        if isinstance(states, Add):
            result = 0            
            for state in states.args:
                result = result + apply_gates(gate**number_of_applications*state, basis)
            states = result
            states = states.expand()

        #if we have a mul, apply gate to each register and multiply result
        elif isinstance(states, Mul):
            #find the Qbits in the Mul
            states = Mul(*states.args)
            for i in range(len(states.args)):
                if isinstance(states.args[i],Qbit):
                    break
            #if we didn't find one, something is wrong
            if not isinstance(states.args[i],Qbit):
                print states
                raise Exception()
 
            #apply the gate the right number of times to this state
            coefficient = Mul(*(states.args[:i]+states.args[i+1:]))
            states = apply_gates(gate**(number_of_applications)*states.args[i], basis)
            states = coefficient*states            
            states = states.expand()
            
        #If we have a single Qbit, apply to this Qbit
        elif isinstance(states, Qbit):
            basis_name = basis.__class__.__name__
            apply_method_name = '_apply_%s' % basis_name  
            apply_method = getattr(gate, apply_method_name)
            states = apply_method(states)
            states = states.expand()
            number_of_applications -= 1
            while number_of_applications > 0:
                states = apply_gates(gate*states)
                number_of_applications -= 1
       
        #if it's not one of those, there is something wrong
        else:
            raise Exception()

    #tack on any coefficients that were there before and simplify
    states = state_coeff*states
    if isinstance(states, (Mul,Add,Pow)):
        states = states.expand()
    return states

#will want to sort stuff in the future so that they are in the correct order                
    
"""
# Look at dimension of basis, only work if it is not a symbol, the dispatch to Qbits._represent_BasisClass
Qbits.represent(self, basis):  
Qbits._represent_XBasisSet(self, dimension):
Qbits._represent_YBasisSet(self, dimension):
"""

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
            
    #if sympy simplified by pulling out a constant coefficeint, undo that
    if isinstance(result, (Mul,Add,Pow)):
        result = result.expand()
    return result

def qbits_to_matrix(qbits):
    #get rid of multiplicative constants
    qbits = qbits.expand()
    
    #if we have a Mul object, find the qbit part qbits to matrix it
    if isinstance(qbits, Mul):
        for i in range(len(qbits.args)):
            if isinstance(qbits.args[i], Qbit):
                break
        if not isinstance(qbits.args[i], Qbit):
            raise Exception()
        #recursively turn qbit into matrix
        return Mul(*(qbits.args[:i] + qbits.args[i+1:]))*qbits_to_matrix(qbits.args[i]) 
    #recursively turn each item in an add into a matrix
    elif isinstance(qbits, Add):
        result = qbits_to_matrix(qbits.args[0])
        for element in qbits.args[1:]:
            result = result + qbits_to_matrix(element)
        return result
    #if we are at the bottom of the recursion, have the base case be representing the matrix
    elif isinstance(qbits, Qbit):
        return qbits._represent()
    else:
        raise Exception("Malformed input")

def represent(circuit, basis = ZBasisSet(), GateRep = False, HilbertSize = None):
    """
        Represents the elements in a certain basis 
    """

    basis_name = basis.__class__.__name__
    rep_method_name = '_represent_%s' % basis_name

    # check if the last element in circuit is Gate
    # if not raise exception becuase size of Hilbert space undefined
    if isinstance(circuit, Qbit):
        return circuit._represent()
    elif isinstance(circuit, Gate):
        if HilbertSize == None:
            raise HilbertSpaceException("User must specify HilbertSize when gates are not applied on Qbits") 
        gate = circuit
        rep_method = getattr(gate, rep_method_name)
        gate_rep = rep_method(HilbertSize)
        return gate_rep
    elif not isinstance(circuit, Mul):
        raise Exception()
    

    qbit = circuit.args[len(circuit.args)-1]
    if isinstance(qbit, Qbit):
        HilbertSize = len(qbit)
        #Turn the definite state of Qbits |X> into a single one in the Xth element of its column vector
        result = qbit._represent()
    elif HilbertSize == None:    
        raise HilbertSpaceException("User must specify HilbertSize when gates are not applied on Qbits")        
    else:
        gate = qbit
        rep_method = getattr(gate, rep_method_name)
        gate_rep = rep_method(HilbertSize)
        result = gate_rep


    #go through each gate (from left->right) and apply it in X basis
    for gate in reversed(circuit.args[:len(circuit.args)-1]):
        basis_name = basis.__class__.__name__
        rep_method_name = '_represent_%s' % basis_name
        if isinstance(gate, Pow):
            number_of_applications = gate.exp
            gate = gate.base
            rep_method = getattr(gate, rep_method_name)
            gate_rep = rep_method(HilbertSize)
            for i in range(number_of_applications):
                result = gate_rep*result
        elif hasattr(gate,  rep_method_name):
            rep_method = getattr(gate, rep_method_name)
            gate_rep = rep_method(HilbertSize)
            result = gate_rep*result
        else:
            raise Exception()     

    return result


def gatesimp(circuit):
    """ will simplify gates symbolically"""
    #Pull gates out of inner Add's and Mul's?

    #bubble sort(?) out gates that commute
    circuit = gatesort(circuit)
    
    #do simplifications
    if isinstance(circuit, Mul):
        for i in range(len(circuit.args)):
            #H,X,Y or Z squared is 1. T**2 = S, S**2 = Z 
            if isinstance(circuit.args[i], Pow):
                if isinstance(circuit.args[i].base, (HadamardGate, XGate, YGate, ZGate)) and isinstance(circuit.args[i].exp, Integer):
                    newargs = (circuit.args[:i] + (circuit.args[i].base**(circuit.args[i].exp % 2),) + circuit.args[i+1:])
                    circuit = gatesimp(Mul(*newargs))
                    break
                elif isinstance(circuit.args[i].base, PhaseGate):
                    newargs = (circuit.args[:i] + (ZGate(circuit.args[i].base.args[0])**(Integer(circuit.args[i].exp/2)), circuit.args[i].base**(circuit.args[i].exp % 2)) + circuit.args[i+1:])
                    circuit =  gatesimp(Mul(*newargs))
                    break
                elif isinstance(circuit.args[i].base,TGate):
                    newargs = (circuit.args[:i] + (SGate(circuit.args[i].base.args[0])**Integer(circuit.args[i].exp/2), circuit.args[i].base**(circuit.args[i].exp % 2)) + circuit.args[i+1:])
                    circuit =  gatesimp(Mul(*newargs))
                    break
            #take care of HYH=-Y,HXH=Z,HZH=X?
            
    return circuit

def gatesort(circuit):
    #recursive bubble sort of gates checking for commutivity
    changes = False
    cirArray = circuit.args
    for i in range(len(cirArray)-1):
        #Go through each element and switch ones that are in wrong order
        if isinstance(cirArray[i], (Gate, Pow)) and isinstance(cirArray[i+1], (Gate, Pow)):
            if isinstance(cirArray[i], Pow):
                first = cirArray[i].base
            else:
                first = cirArray[i]

            if isinstance(cirArray[i+1], Pow):
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
                    circuit = Mul(*(circuit.args[:i] + (circuit.args[i+1],) + (circuit.args[i],) + circuit.args[i+2:])) 
                    cirArray = circuit.args
                    changes = True
    #if we made changes, recursively call the sort method until we don't
    if changes:
        return gatesort(circuit)
    else:
        return circuit    

class HilbertSpaceException(Exception):
    pass

