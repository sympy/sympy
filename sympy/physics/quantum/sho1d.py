"""Simple Harmonic Oscillator 1-Dimension"""

from sympy import sqrt, I, Symbol, Integer
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.state import Bra, Ket, State
from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.cartesian import X, Px
from sympy.functions.special.tensor_functions import KroneckerDelta

#-------------------------------------------------------------------------

class SHOOp(Operator):
	"""A base class for the SHO Operators.
	
	We are limiting the number of arguments to be 1.
	
	"""
	
	@classmethod
	def _eval_args(cls, args):
		args = QExpr._eval_args(args)
		if len(args) == 1:
			return args
		else:
			raise ValueError("Too many arguments")

			
class RaisingOp(SHOOp):
	"""The Raising Operator or a^dagger.
	
	When a^dagger acts on a state it raises the state up by one. Taking 
	the adjoint of a^dagger returns 'a', the Lowering Operator. a^dagger 
	be rewriten in terms of postion and momentum.
	
	Parameters
	==========
	
	args : tuple
		The list of numbers or parameters that uniquely specify the 
		operator.
	
	Examples
	========
	
	Create a Raising Operator and rewrite it in terms of positon and
	momentum, and show that taking its adjoint returns 'a':
	
		>>> from sympy.physics.quantum import sho1d
		>>> from sympy.functions.elementary.complexes import adjoint
	
		>>> ad = RaisingOp('a')
		>>> ad().rewrite('xp').doit()
		1/sqrt(2*hbar*m*omega)*(-I*Px + m*omega*X)
		
		>>> adjoint(ad)
		a
		
	Taking the commutator of a^dagger with other Operators:
	
		>>> from sympy.physics.quantum import sho1d, Commutator
		
		>>> ad = RaisingOp('a')
		>>> a = LoweringOp('a')
		>>> N = NumberOp('N')
		>>> Commutator(ad, a).doit()
		-1
		>>> Commutator(ad, N).doit()
		-a^dagger
	
	Apply a^dagger to a state:
	
		>>> from sympy.physics.quantum import sho1d, qapply
		
		>>> ad = RaisingOp('a')
		>>> k = SHOKet('k')
		>>> qapply(ad,k)
		sqrt(k+1)|k+1>
		
	"""

	def _eval_rewrite_as_xp(self, *args):
		return (Integer(1)/sqrt(Integer(2)*hbar*m*omega))*(
			Integer(-1)*I*Px + m*omega*X)

	def _eval_adjoint(self):
		return LoweringOp(*self.args)
		
	def _eval_commutator_LoweringOp(self, other):
		return Integer(-1)
	
	def _eval_commutator_NumberOp(self, other):
		return Integer(-1)*self
		
	def _apply_operator_SHOKet(self, ket):
		temp = ket.n + Integer(1)
		return sqrt(temp)*SHOKet(temp)
	
	#---------------------------------------------------------------------
	# Printing Methods
	#---------------------------------------------------------------------
	
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
	"""The Lowering Operator or 'a'.
	
	When 'a' acts on a state it lowers the state up by one. Taking 
	the adjoint of 'a' returns a^dagger, the Raising Operator. 'a' 
	be rewriten in terms of postion and momentum.
	
	Parameters
	==========
	
	args : tuple
		The list of numbers or parameters that uniquely specify the 
		operator.
	
	Examples
	========
	
	Create a Lowering Operator and rewrite it in terms of positon and
	momentum, and show that taking its adjoint returns a^dagger:
	
		>>> from sympy.physics.quantum import sho1d
		>>> from sympy.functions.elementary.complexes import adjoint
	
		>>> a = LoweringOp('a')
		>>> a().rewrite('xp').doit()
		1/sqrt(2*hbar*m*omega)*(I*Px + m*omega*X)
		
		>>> adjoint(a)
		a^dagger
		
	Taking the commutator of 'a' with other Operators:
	
		>>> from sympy.physics.quantum import sho1d, Commutator
		
		>>> a = LoweringOp('a')
		>>> ad = RaisingOp('a')
		>>> N = NumberOp('N')
		>>> Commutator(a, ad).doit()
		1
		>>> Commutator(a, N).doit()
		a
	
	Apply 'a' to a state:
	
		>>> from sympy.physics.quantum import sho1d, qapply
		
		>>> a = LoweringOp('a')
		>>> k = SHOKet('k')
		>>> qapply(a,k)
		sqrt(k)|k-1>
		
	Taking 'a' of the lowest state will return 0:
	
		>>> from sympy.physics.quantum import sho1d, qapply
		
		>>> a = LoweringOp('a')
		>>> k = SHOKet(0)
		>>> qapply(a,k)
		0
		
	"""
	
	def _eval_rewrite_as_xp(self, *args):
		return (Integer(1)/sqrt(Integer(2)*hbar*m*omega))*(
			I*Px + m*omega*X)
	
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
	"""The Number Operator is simply a^dagger*a
	
	It is often useful to write a^dagger*a as simply the Number Operator
	because the Number Operator commutes with the Hamiltonian. And can be
	expressed using the Number Operator. Also the Number Operator can be
	applied to states.
	
	Parameters
	==========
	
	args : tuple
		The list of numbers or parameters that uniquely specify the 
		operator.
	
	Examples
	========
	
	Create a Number Operator and rewrite it in terms of the ladder 
	operators and Hamiltonian:
	
		>>> from sympy.physics.quantum import sho1d
		
		>>> N = NumberOp('N')
		>>> N().rewrite('a').doit()
		a^dagger*a
		>>> N().rewrite('H').doit()
		H/(hbar*omega) - 1/2
		
	Take the Commutator of the Number Operator with other Operators:
	
		>>> from sympy.physics.quantum import sho1d, Commutator
		
		>>> N = NumberOp('N')
		>>> H = Hamiltonian('H')
		>>> ad = RaisingOp('a')
		>>> a = LoweringOp('a')
		>>> Commutator(N,H).doit()
		0
		>>> Commutator(N,ad).doit(0
		a^dagger
		>>> Commutator(N,a).doit()
		-a
		
	Apply the Number Operator to a state:
	
		>>> from sympy.physics.quantum import sho1d, qapply
		
		>>> N = NumberOp('N')
		>>> k = SHOKet('k')
		>>> qapply(N,k)
		k|k>
	
	"""
	
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
	"""The Hamiltonian Operator.
	
	The Hamiltonian is used to solve the time-independent Schrodinger
	equation. The Hamiltonian can be expressed using the ladder operators,
	as well as by position and momentum.
	
	Parameters
	==========
	
	args : tuple
		The list of numbers or parameters that uniquely specify the 
		operator.
	
	Examples
	========
	
	Create a Hamiltonian Operator and rewrite it in terms of the ladder
	operators, position and momentum, and the Number Operator:
	
		>>> from sympy.physics.quantum import sho1d
		
		>>> H = Hamiltonian('H')
		>>> H().rewrite('a').doit()
		hbar*omega*(a^dagger*a + 1/2)
		>>> H().rewrite('xp').doit()
		1/2m * (Px**2 + (m*omega*X)^2)
		>>> H().rewrite('n').doit()
		hbar*omega*(N + 1/2)
		
	Take the Commutator of the Hamiltonian and the Number Operator:
	
		>>> from sympy.physics.quantum import sho1d, Commutator
		
		>>> H = Hamiltonian('H')
		>>> N = NumberOp('N')
		>>> Commutator(H,N).doit()
		0
		
	Apply the Hamiltonian Operator to a state:
	
		>>> from sympy.physics.quantum import sho1d, qapply
		
		>>> H = Hamiltonian('H')
		>>> k = SHOKet('k')
		>>> qapply(H,k)
		hbar*omega*(k + 1/2)|k>
	
	"""

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
	
#-------------------------------------------------------------------------

class SHOState(State):
	"""State class for SHO states"""
	
	@property
	def n(self):
		return self.args[0]


class SHOKet(SHOState, Ket):
	"""1D eigenket.
	
	Inherits from SHOState and Ket.
	
	Parameters
	==========
	
	args : tuple
		The list of numbers or parameters that uniquely specify the ket.
		This is usually its quantum numbers or its symbol.
	
	Examples
	========
	
	Ket's know about their associated bra:
	
		>>> from sympy.physics.quantum import sho1d
		
		>>> k = SHOKet('k')
		>>> k.dual
		<k|
		>>> k.dual_class()
		<class 'sympy.quantum.sho1d.SHOBra'>
		
	Take the Inner Product with a bra:
	
		>>> from sympy.physics.quantum import sho1d, InnerProduct
		
		>>> k = SHOKet('k')
		>>> b = SHOBra('b')
		>>> InnerProduct(k,b).doit()
		<b|k>
		
	"""
	
	@classmethod
	def dual_class(self):
		return SHOBra

	def _eval_innerproduct_SHOBra(self, bra, **hints):
		result = KroneckerDelta(self.n, bra.n)
		return result
	
		
class SHOBra(SHOState, Bra):
	"""A time-independent Bra in SHO.
	
	Inherits from SHOState and Bra.
	
	Parameters
	==========
	
	args : tuple
		The list of numbers or parameters that uniquely specify the ket.
		This is usually its quantum numbers or its symbol.
	
	Examples
	========
	
	Bra's know about their associated ket:
	
		>>> from sympy.physics.quantum import sho1d
		
		>>> b = SHOBra('b')
		>>> b.dual
		|b>
		>>> b.dual_class()
		<class 'sympy.quantum.sho1d.SHOKet'>
		
	"""
	
	@classmethod
	def dual_class(self):
		return SHOKet

		
ad = RaisingOp('a')
a = LoweringOp('a')
H = Hamiltonian('H')
N = NumberOp('N')
omega = Symbol('omega')
m = Symbol('m')