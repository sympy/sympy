"""Simple Harmonic Oscillator 1-Dimension"""

#%load_ext sympyprinting

from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.qexpr import *
from sympy.physics.quantum.cartesian import *

#--------------------------------------------------------------------

class SHOOp(Operator):
	
	@classmethod
	def _eval_args(cls, args):
		args = QExpr._eval_args(args)
		if len(args) == 1:
			return args
		else:
			raise ValueError("Too many arguments")

			
class RaisingOp(SHOOp):

	def _eval_rewrite_as_xp(self, *args):
		return (Integer(1)/sqrt(Integer(2)*hbar*m*w))*(Integer(-1)*I*Px + m*w*X)

	def _eval_adjoint(self):
		return LoweringOp(*self.args)
		
	def _eval_commutator_LoweringOp(self, other):
		return Integer(-1)
	
	def _eval_commutator_NumberOp(self, other):
		return Integer(-1)*self
		
	def _apply_operator_SHOKet(self, ket):
		temp = ket.n + Integer(1)
		return sqrt(temp)*SHOKet(temp)
		
	# Printing Methods
	def _sympyrepr(self, printer, *args):
		arg0 = printer._print(self.args[0], *args)
		return '%s(%s)' % (self.__class__.__name__, arg0)
		
	def _sympystr(self, printer, *args):
		arg0 = printer._print(self.args[0], *args)
		return '%s(%s)' % (self.__class__.__name__, arg0)
		
	def _pretty(self, printer, *args):
		from sympy.printing.pretty.stringpict import prettyForm
		pform = printer._print(self.args[0], *args)
		pform = pform**prettyForm(u'\u2020')
		return pform
		
	def _latex(self, printer, *args):
		arg = printer._print(self.args[0])
		return '%s^{\\dag}' % arg
	
	
	def _print_label_repr(self, printer, *args):
		return self._print_sequence(
			self.label, self._sympyrepr, ',',printer, *args
		)
		
	def _print_label_pretty(self,printer, *args):
		return self._print_sequence(
			self.label, self._pretty, self._label_separator, printer, *args
		)
		
	def _print_label_latex(self, printer, *args):
		return self._print_sequence(
			self.label, self._latex, self._label_separator, printer, *args
		)

		
class LoweringOp(SHOOp):

	def _eval_rewrite_as_xp(self, *args):
		return (Integer(1)/sqrt(Integer(2)*hbar*m*w))*(I*Px + m*w*X)
	
	def _eval_adjoint(self):
		return RaisingOp(*self.args)

	def _eval_commutator_RaisingOp(self, other):
		return Integer(1)
		
	def _eval_commutator_NumberOp(self, other):
		return Integer(1)*self

	def _apply_operator_SHOKet(self, ket):
		temp = ket.n - Integer(1)
		if temp == Integer(0):
			return Integer(0)
		else:
			return sqrt(ket.n)*SHOKet(temp)


class NumberOp(SHOOp):
	
	def _eval_rewrite_as_a(self, *args):
		return ap*am
	
	def _eval_rewrite_as_H(self, *args):
		return H/(hbar*w) - Integer(1)/Integer(2)
	
	@classmethod
	def _apply_operator_SHOKet(self, ket):
		return ket.n*ket
		
	def _eval_commutator_Hamiltonian(self, other):
		return Integer(0)
		
	def _eval_commutator_RaisingOp(self, other):
		return other
		
	def _eval_commutator_LoweringOp(self, other):
		return Integer(-1)*other


class Hamiltonian(SHOOp):

	def _eval_rewrite_as_a(self, *args):
		return hbar*w*(ap*am + Integer(1)/Integer(2))
		
	def _eval_rewrite_as_xp(self, *args):
		return (Integer(1)/(Integer(2)*m))*(Px**2 + (m*w*X)**2)
		
	def _eval_rewrite_as_n(self, *args):
		return hbar*w*(N + Integer(1)/Integer(2))
		
	def _eval_rewrite_as_am(self, *args):
		return hbar*w*(am*ap - Integer(1)/Integer(2))
	
	@classmethod
	def _apply_operator_SHOKet(self, ket):
		return (hbar*w*(ket.n + Integer(1)/Integer(2)))*ket
		
	def _eval_commutator_NumberOp(self, other):
		return Integer(0)
	
#--------------------------------------------------------------------

class SHOState(State):

	@property
	def n(self):
		return self.args[0]


class SHOKet(SHOState, Ket):
	
	@classmethod
	def dual_class(self):
		return SHOBra

	def _eval_innerproduct_SHOBra(self, bra, **hints):
		result = KroneckerDelta(self.n, bra.n)
		return result
	
		
class SHOBra(SHOState, Bra):

	@classmethod
	def dual_class(self):
		return SHOKet

		
ap = RaisingOp('a')
am = LoweringOp('a')
H = Hamiltonian('H')
N = NumberOp('N')
k = SHOKet('k')
b = SHOBra('b')
w = Symbol('omega')
m = Symbol('m')