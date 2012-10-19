"""Simple Harmonic Oscillator 1-Dimension"""

%load_ext sympyprinting

from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.qexpr import *
#from sympy.physics.quantum.cartesian import*
#from sympy.printing.pretty.stringpict import prettyForm

#--------------------------------------------------------------------

class RaisingOp(Operator):
	
	@classmethod
	def default_args(self):
		return("a",)
	
	def _eval_adjoint(self):
		return Dagger(self)
		
	# Printing
	
	_label_separator = ','
	
	def _print_sequence(self, seq, printer, *args):
		result = []
		for item in seq:
			result.append(printer._print(item, *args))
		return sep.join(result)
		
	def _print_sequence_pretty(self, seq, printer, *args):
		pform = printer._print(seq[0], *args)
		for item in seq[1]:
			pform = prettyForm(*pform.right((sep)))
			pform = prettyForm(*pform.right((printer._print(item, *args))))
		return pform
	
	def _print_label(self, printer, *args):
		return self._print_sequence(
			self.label, self._label_separator, printer, *args
		)
	
	def _print_label_repr(self, printer, *args):
		return self._print_sequence(
			self.label, ',', printer, *args
		)
		
	def _print_label_pretty(self, printer, *args):
		return self._print_sequence_pretty(
			self.label, self._label_separator, printer, *args
		)
		
	def _print_label_latex(self, printer, *args):
		return self._print_sequence(
			self.label, self._label_separator, printer, *args
		)
		
	def _print_contents(self, printer, *args):
		return self._printer_label(printer, *args)
		
	def _print_contents_repr(self, printer, *args):
		return self._print_label_repr(printer, *args)
		
	def _print_contents_pretty(self, printer, *args):
		return self._print_label_pretty(printer, *args)
		
	def _print_contents_latex(self, printer, *args):
		return self._print_label_latex(printer, *args)

class LoweringOp(Operator):

	@classmethod
	def default_args(self):
		return ("a",)
		
	def _eval_adjoint(self):
		return Dagger(self)


#--------------------------------------------------------------------

		
class SHOBase(StateBase):
	pass


	
class SHOState(SHOBase):
	pass


class SHOKet(SHOState, KetBase):
	
	@classmethod
	def dual_class(self):
		return SHOBra
		
class SHOBra(SHOState, BraBase):

	@classmethod
	def dual_class(self):
		return SHOKet