from sympy.physics.quantum import Operator, State, Ket, Bra
from sympy.physics.hilbert import HilbertSpace
from sympy.physics.units import hbar, kg, m, planck
from sympy.core.numbers import Pi
from sympy import Rational, Symbol, I
from sympy.printing.str import sstr
from sympy.functions.elementary.trigonometric import sin
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.exponential import exp

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

    def time_dep_eigenvector(self, n, t=Symbol('t')):
        return TimeDepOpPIAB('U').rep(n, t)*PIABKet(self.name, n)

    @property
    def eigenvalue(self, n):
        return self.n**2*Pi()**2*hbar**2/(2*self.mass*(self.L)**2)

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

class TimeDepOpPIAB(PIAB, Operator):

    def rep(self, n, t=Symbol('t')):
        return exp(I*t*self.energy(n)/planck)

    def energy(self, n):
        return n**2*Pi()**2*hbar**2/(2*self.mass*(self.L)**2)
