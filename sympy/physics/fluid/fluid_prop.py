import numpy as np
from sympy.core import S, Symbol, diff

	'''
	Density: The density, or more precisely, the volumetric mass density, of a substance is its mass per unit volume
	================================================================================================================
	Input Paramters :
	 m - Mass of the Fluid
	 v - Volume of the Fluid

	Output Parameter :
	 d - Density of the Fluid

	================================================================================================================
	Example :

	>>> from sympy.core import S, Symbol, symbols
	>>> from sympy.physics.fluid import *
	>>> m, v = symbols('m, v')
	>>> d = density(m, v)
	'''

def density(mass, volume):
	return mass/volume;
