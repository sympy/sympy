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
		return LoweringOp(*self.args)
		
	# Printing look at Dagger code
	
class LoweringOp(Operator):

	@classmethod
	def default_args(self):
		return ("a",)
		
	def _eval_adjoint(self):
		return RaisingOp(*self.args)


#--------------------------------------------------------------------

		
	
class SHOState(State):
	pass


class SHOKet(SHOState, Ket):
	
	@classmethod
	def dual_class(self):
		return SHOBra
		
class SHOBra(SHOState, Bra):

	@classmethod
	def dual_class(self):
		return SHOKet
		
		
#Make ipython notebooks to test