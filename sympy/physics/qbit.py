"""
    Single qbits and their gates
"""
from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify, Matrix, elementary
from sympy.core.basic import S, sympify
from sympy.core.function import Function

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
    An object representing a Hadamard Gate:
    1/sqrt(2)*[[1, 1], [1, -1]]
    """
    def _sympystr(self, printer, *args):
        return "H(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1., 1.], [1., -1.]])*(1./(2.**(1./2.)))

    def _represent_XBasisSet(self, HilbertSize):
        pass

    def _represent_YBasisSet(self, HilbertSize):
        pass

class XGate(Gate):
    """
    An object representing a Pauli-X gate:
    [[0, 1], [1, 0]]
    """
    def _sympystr(self, printer, *args):
        return "X(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[0, 1.], [1., 0]])

    def _represent_XBasisSet(self, HilbertSize):
        pass

    def _represent_YBasisSet(self, HilbertSize):
        pass

class YGate(Gate):
    """
    An object representing a Pauli-Y gate:
    [[0, -i], [i, 0]]
    """
    def _sympystr(self, printer, *args):
        return "Y(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[0, -complex(0,1)], [complex(0,1), -0]])

    def _represent_XBasisSet(self, HilbertSize):
        pass

    def _represent_YBasisSet(self, HilbertSize):
        pass

class ZGate(Gate):
    """
    An object representing a Pauli-Z gate:
    [[1, 0], [0, -1]]
    """
    def _sympystr(self, printer, *args):
        return "Z(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1., 0], [0, -1.]])

    def _represent_XBasisSet(self, HilbertSize):
        pass

    def _represent_YBasisSet(self, HilbertSize):
        pass

class PhaseGate(Gate):
    """
    An object representing a phase gate:
    [[1, 0], [0, i]]
    """
    def _sympystr(self, printer, *args):
        return "S(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1, 0], [0, i]])

    def _represent_XBasisSet(self, HilbertSize):
        pass

    def _represent_YBasisSet(self, HilbertSize):
        pass

class TGate(Gate):
    """
    An object representing a pi/8 gate:
    [[1, 0], [0, e**(i*pi/4)]]
    """
    def _sympystr(self, printer, *args):
        return "T(%s)" % printer._print(self.args[0], *args)

    def _represent_ZBasisSet(self, HilbertSize):
        return Matrix([[1., 0], [1., -1.]])

    def _represent_XBasisSet(self, HilbertSize):
        pass

    def _represent_YBasisSet(self, HilbertSize):
        pass

def represent(circuit, basis, GateRep = False):
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

    #for now, assume Qbits are all zero at beginning (just for testing).
    result = Matrix([1, 0])


    #go through each gate (from left->right) and apply it in X basis
    #http://docs.sympy.org/modules/mpmath/matrices.html says:
    #"Matrices in mpmath are implemented using dictionaries. Only non-zero values are stored, so it is cheap to represent sparse matrices."
    for gate in reversed(circuit.args[:len(circuit.args)-1]):
        basis_name = basis.__class__.__name__
        rep_method_name = '_represent_%s' % basis_name 
        if hasattr(gate,  rep_method_name):
            rep_method = getattr(gate, rep_method_name)
            gate_rep = rep_method(HilbertSize)
            result = gate_rep*result        

    print result

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

