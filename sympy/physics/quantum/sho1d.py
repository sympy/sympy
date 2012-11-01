"""Simple Harmonic Oscillator 1-Dimension"""

from sympy import sqrt, I, Symbol, Integer
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.state import Bra, Ket, State
from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.cartesian import X, Px
from sympy.functions.special.tensor_functions import KroneckerDelta

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
		return (Integer(1)/sqrt(Integer(2)*hbar*m*omega))*(Integer(-1)*I*Px + m*omega*X)

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
		return (Integer(1)/sqrt(Integer(2)*hbar*m*omega))*(I*Px + m*omega*X)
	
	def _eval_adjoint(self):
		return RaisingOp(*self.args)

	def _eval_commutator_RaisingOp(self, other):
		return Integer(1)
		
	def _eval_commutator_NumberOp(self, other):
		return Integer(1)*self

	def _apply_operator_SHOKet(self, ket):
		temp = ket.n - Integer(1)
		if ket.n == Integer(0):
			return Integer(0)
		else:
			return sqrt(ket.n)*SHOKet(temp)


class NumberOp(SHOOp):
	
	def _eval_rewrite_as_a(self, *args):
		return ad*a
	
	def _eval_rewrite_as_H(self, *args):
		return H/(hbar*omega) - Integer(1)/Integer(2)
	
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
		return hbar*omega*(ad*a + Integer(1)/Integer(2))
		
	def _eval_rewrite_as_xp(self, *args):
		return (Integer(1)/(Integer(2)*m))*(Px**2 + (m*omega*X)**2)
		
	def _eval_rewrite_as_n(self, *args):
		return hbar*omega*(N + Integer(1)/Integer(2))
	
	def _apply_operator_SHOKet(self, ket):
		return (hbar*omega*(ket.n + Integer(1)/Integer(2)))*ket
		
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

		
ad = RaisingOp('a')
a = LoweringOp('a')
H = Hamiltonian('H')
N = NumberOp('N')
omega = Symbol('omega')
m = Symbol('m')