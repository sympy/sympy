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

'''
	Specific Weight: The specific weight (also known as the unit weight) is the weight per unit volume of a material.
	================================================================================================================
	Input Paramters :
	 m - Mass of the Fluid
	 v - Volume of the Fluid
	 g - Acceleration due to Gravity

	Output Parameter :
	 w - Specific Weight

	================================================================================================================
	Example :

	>>> from sympy.core import S, Symbol, symbols
	>>> from sympy.physics.fluid import *
	>>> m, v, g = symbols('m, v, g')
	>>> w = sp_weight(m, v, g)
'''

def sp_weight(mass, volume , gravity):
	return (mass*gravity)/volume;
