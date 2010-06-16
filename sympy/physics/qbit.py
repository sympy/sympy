"""
    Single qbits and their gates
"""
from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify


class Qbit(Expr):
    """
    Represents a single quantum gate
    """
    def __new__(cls, *args):
        obj = Expr.__new__(cls, args, commutative = False)
        return obj

    def dimension(self):
        return len(self.args[0])

    def _sympyrepr_(self):
        pass

class Gate(Expr):
    """
    A gate operator that acts qubit(s)
    """
    def __new__(cls, *args):
        obj = Expr.__new__(cls, args, commutative = False)
        return obj

    def _sympyrepr_(self, printer, *args):
        return "%s(%i)" %  (self.__class__.__name__, printer._print(self.name, *args))


class HadamardGate(Gate):
    """
    An object representing a Hadamard Gate:
    1/sqrt(2)*[[1, 1], [1, -1]]
    """
    def _sympystr_(self, printer, *args):
        return printer._print("H(%i)", *args)
   
class XGate(Gate):
    """
    An object representing a Pauli-X gate:
    [[0, 1], [1, 0]]
    """
    def __str__(self):
        return "X(" + str(self.args[0][0]) + ")" 
    
class YGate(Gate):
    """
    An object representing a Pauli-Y gate:
    [[0, -i], [i, 0]]
    """
    def __str__(self):
        return "Y(" + str(self.args[0][0]) + ")" 
    
class ZGate(Gate):
    """
    An object representing a Pauli-Z gate:
    [[1, 0], [0, -1]]
    """
    def __str__(self):
        return "Z(" + str(self.args[0][0]) + ")" 

class PhaseGate(Gate):
    """
    An object representing a phase gate:
    [[1, 0], [0, i]]
    """
    def __str__(self):
        return "S(" + str(self.args[0][0]) + ")" 

class TGate(Gate):
    """
    An object representing a pi/8 gate:
    [[1, 0], [0, e**(i*pi/4)]]
    """ 
    def __str__(self):
        return "T(" + str(self.args[0][0]) + ")" 
