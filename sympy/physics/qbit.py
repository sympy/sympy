"""
    Single qbits and their gates
    a = [[x for x in 2*q] for q in l]+[[x for x in 2*q] for q in l]
"""
from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify, Matrix, elementary
from sympy.core.numbers import *
from sympy.core.basic import S, sympify
from sympy.core.function import Function
from sympy.functions.elementary.exponential import *
from sympy.functions.elementary.miscellaneous import *

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

    def _sympystr(self, printer, *args):
        string = ""
        for it in self.args:
            string = string + str(it)
        return "|%s>" % printer._print(string, *args)

    def _sympyrepr(self, printer, *args):
        return "%s(%i)" %  (self.__class__.__name__, printer._print(self.name, args))

class QbitSpace(Expr):
    def __new__(cls, size):
        return Expr.__new__(cls, size, commutative = False)

    def _sympystr(self,printer, *args):
        return "l2(2)**%d" % (self.args[0])

    @property
    def dimension(self):
        return size

class Gate(Expr):
    """
    A gate operator that acts on qubit(s)
    (will need tensor product to get it working for multiple Qbits)
    """
    def __new__(cls, *args):
        obj = Expr.__new__(cls, *args, commutative = False)
        return obj

    @property
    def minimumdimension(self):
        return max(self.args)

    def _sympyrepr(self, printer, *args):
        return "%s(%s)" %  (self.__class__.__name__, printer._print(self.name, args))

    def _represent_ZBasisSet(self, HilbertSize):
        raise NotImplementedError("Z-Basis Representation not implemented")

    def _represent_XBasisSet(self, HilbertSize):
        raise NotImplementedError("X-Basis Representation not implemented")

    def _represent_YBasisSet(self, HilbertSize):
        raise NotImplementedError("Y-Basis Representation not implemented")

def representHilbertSpace(gate, HilbertSize, qbit):
    I = Matrix([[1,0],[0,1]])
    product = [I for x in range(HilbertSize)]
    product[qbit] = gate
    MatrixRep = TensorProduct(product)
    return MatrixRep

def TensorProduct(args):
    MatrixExpansion = args[len(args)-1]
    for gate in args[:len(args)-1]:
        a = ((gate[0]*MatrixExpansion).tolist())
        b = ((gate[1]*MatrixExpansion).tolist())
        c = ((gate[2]*MatrixExpansion).tolist())
        d = ((gate[3]*MatrixExpansion).tolist())
        b = [a[x].extend(b[x]) for x in range(len(a))]
        d = [c[x].extend(d[x]) for x in range(len(c))]
        MatrixExpansion = Matrix(a + c)
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

    def _represent_ZBasisSet(self, HilbertSize):
        if self.minimumdimension > HilbertSize:
            raise HilbertSpaceException()
        gate = Matrix([[1, 1], [1, -1]])*(1/sqrt(2))
        if HilbertSize  == 1:            
            return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))
        else:
            m = representHilbertSpace(gate, HilbertSize, self.args[0])
            return m

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[1, 1], [1, -1]])*(1/sqrt(2))

    def _represent_YBasisSet(self, HilbertSize):
        pass

class XGate(Gate):
    """
    An object representing a Pauli-X gate:
    """
    def _sympystr(self, printer, *args):
        return "X(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        if self.minimumdimension > HilbertSize:
            raise HilbertSpaceException()
        return Matrix([[0, 1], [1, 0]])

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[1,0],[0,-1]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

class YGate(Gate):
    """
    An object representing a Pauli-Y gate:
    """
    def _sympystr(self, printer, *args):
        return "Y(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[0, -complex(0,1)], [complex(0,1), -0]])

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[0,complex(0,1)],[complex(0,-1),0]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

class ZGate(Gate):
    """
    An object representing a Pauli-Z gate:
    """
    def _sympystr(self, printer, *args):
        return "Z(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1, 0], [0, -1]])

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[0,1],[1,0]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

class PhaseGate(Gate):
    """
    An object representing a phase gate:
    """
    def _sympystr(self, printer, *args):
        return "S(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1, 0], [0, complex(0,1)]])

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[complex(.5,.5), complex(.5,-.5)], [complex(.5,-.5),complex(.5,.5)]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

class TGate(Gate):
    """
    An object representing a pi/8 gate:
    """
    def _sympystr(self, printer, *args):
        return "T(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1, 0], [0, exp(complex(0,Pi/4))]])

    def _represent_XBasisSet(self, HilbertSize):
        return Matrix([[.5+.5*exp(complex(0,Pi/4)),.5-.5*exp(complex(0,Pi/4))],[.5-.5*exp(complex(0,Pi/4)),.5+.5*exp(complex(0,Pi/4))]])

    def _represent_YBasisSet(self, HilbertSize):
        pass

def represent(circuit, basis = ZBasisSet(), GateRep = False):
    """
        Represents the elements in a certain basis 
    """

    # check if the last element in circuit is Gate
    # if not raise exception becuase size of Hilbert space undefined
    wires = circuit.args[len(circuit.args)-1]    
    if isinstance(wires, Qbit):
        HilbertSize = len(wires)
    else:
        raise HilbertSpaceException()

    #Turn the definite state of Qbits |X> into a single one in the Xth element of its column vector
    n = 1
    definiteState = 0
    for it in wires.args:
        definiteState += n*it
        n = n*2
    result = [0 for x in range(2**HilbertSize)]
    result[definiteState] = 1

    result = Matrix(result)

    #go through each gate (from left->right) and apply it in X basis
    #http://docs.sympy.org/modules/mpmath/matrices.html says:
    #"Matrices in mpmath are implemented using dictionaries. Only non-zero values are stored, so it is cheap to represent sparse matrices."
    for gate in reversed(circuit.args[:len(circuit.args)-1]):
        basis_name = basis.__class__.__name__
        rep_method_name = '_represent_%s' % basis_name 
        if hasattr(gate,  rep_method_name):
            rep_method = getattr(gate, rep_method_name)
            gate_rep = rep_method(HilbertSize)
            print gate_rep
            result = gate_rep*result        

    return result

"""
#find the basis it should be represented by for each gate
for gate in gates:
    basis_name = basis.__class__.__name__
    rep_method_name = '_represent_%s' % basis_name 
    if hasattr(gate,  rep_method_name):
        rep_method = getattr(gate, rep_method_name)
        gate_rep = rep_method(size)
"""

class HilbertSpaceException(Exception):
    pass

