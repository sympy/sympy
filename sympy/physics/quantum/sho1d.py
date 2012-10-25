"""Simple Harmonic Oscillator 1-Dimension"""


#%load_ext sympyprinting

from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.qexpr import *
#from sympy.physics.quantum.cartesian import*

#--------------------------------------------------------------------

class RaisingOp(Operator):
	
	@classmethod
	def default_args(self):
		args = QExpr._eval_args(args)
		if len(args) == 1:
			return("a",)
		else:
			raise ValueError("Too many arguments")
	
	def _eval_adjoint(self):
		return LoweringOp(*self.args)
		
	def _eval_commutator_LoweringOp(self, other):
		return -1
		
	# Printing 
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

class LoweringOp(Operator):

	@classmethod
	def default_args(self):

		args = QExpr._eval_args(args)
		if len(args) == 1:
			return("a",)
		else:
			raise ValueError("Too many arguments")
		
	def _eval_adjoint(self):
		return RaisingOp(*self.args)
		
	def _eval_commutator_RaisingOp(self, other):
		return 1



#--------------------------------------------------------------------


class SHOState(SHOBase):
	pass


class SHOKet(SHOState, Ket):
	
	@classmethod
	def dual_class(self):
		return SHOBra

	def _eval_innerproduct_SHOBra(self, bra, **hints):
		result = KroneckerDelta(self, bra)
		return result
		
class SHOBra(SHOState, Bra):

	@classmethod
	def dual_class(self):
		return SHOKet

