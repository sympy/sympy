"""Simple Harmonic Oscillator 1-Dimension"""

<<<<<<< HEAD
#%load_ext sympyprinting
=======
%load_ext sympyprinting
>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e

from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.qexpr import *
#from sympy.physics.quantum.cartesian import*
<<<<<<< HEAD

=======
#from sympy.printing.pretty.stringpict import prettyForm
>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e

#--------------------------------------------------------------------

class RaisingOp(Operator):
	
	@classmethod
	def default_args(self):
<<<<<<< HEAD
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
=======
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
>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e
		)
		
	def _print_label_latex(self, printer, *args):
		return self._print_sequence(
<<<<<<< HEAD
			self.label, self._latex, self._label_separator, printer, *args
		)
	
	
	
=======
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

>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e
class LoweringOp(Operator):

	@classmethod
	def default_args(self):
<<<<<<< HEAD
		args = QExpr._eval_args(args)
		if len(args) == 1:
			return("a",)
		else:
			raise ValueError("Too many arguments")
		
	def _eval_adjoint(self):
		return RaisingOp(*self.args)
		
	def _eval_commutator_RaisingOp(self, other):
		return 1
=======
		return ("a",)
		
	def _eval_adjoint(self):
		return Dagger(self)
>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e


#--------------------------------------------------------------------

		
<<<<<<< HEAD
	
class SHOState(State):
	pass


class SHOKet(SHOState, Ket):
=======
class SHOBase(StateBase):
	pass


	
class SHOState(SHOBase):
	pass


class SHOKet(SHOState, KetBase):
>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e
	
	@classmethod
	def dual_class(self):
		return SHOBra
		
<<<<<<< HEAD
	def _eval_innerproduct_SHOBra(self, bra, **hints):
		result = KroneckerDelta(self, bra)
		return result
		
class SHOBra(SHOState, Bra):

	@classmethod
	def dual_class(self):
		return SHOKet
	
=======
class SHOBra(SHOState, BraBase):

	@classmethod
	def dual_class(self):
		return SHOKet
>>>>>>> a9a035523a4ce79aa10c68494d2655b7ef5ae39e
