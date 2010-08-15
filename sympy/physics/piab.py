from sympy.physics.quantum import Operator, State, Ket, Bra
from sympy.physics.hilbert import HilbertSpace
from sympy.physics.units import hbar, kg, m
from sympy.core.numbers import Pi
from sympy import Rational, Symbol
from sympy.printing.str import sstr
from sympy.functions.elementary.trigonometric import sin
from sympy.functions.elementary.miscellaneous import sqrt

me = mass_electron = Rational('9.10938215')*10**-31 * kg #mass of electron

class RHS1D(HilbertSpace):
    pass

class PIAB(object):
    mass = me
    L = 1 * 10**-9 * m

class PIABHamiltonian(Operator, PIAB):
    def __new__(cls, name):
        return Operator.__new__(cls, name)

    def eigenvector(self, n):
        return PIABKet(self.name, n)

class PIABState(State, PIAB):
    __slots__ = ['n']
    
    def __new__(cls, name, n):
        obj = State.__new__(cls, name)
        obj.n = n
        return obj

    @property
    def name(self):        
        return self.args[0]
        
    @property
    def Hamiltonian(self):
        return PIABHamiltonian()

    @property
    def eigenvalue(self):
        return self.n**2*Pi()**2*hbar**2/(2*self.mass*(self.L)**2)

    @property
    def represent_positionbasis(self):
        x = Symbol('x')
        return sqrt(2/(self.L))*sin(self.n*Pi()*x*m/self.L)

class PIABKet(PIABState, Ket):

    @property
    def dual(self):
        return PIABBra(*self.args) 

class PIABBra(PIABState, Bra):

    @property
    def dual(self):
        return PIABKet(*self.args)
