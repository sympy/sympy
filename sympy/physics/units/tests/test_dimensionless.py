from sympy.physics.units import Dimension
from sympy.physics.units import DimensionSystem

""" 
	We need to define the base quantities globally to work with 
	and define the derived quantities from there on
"""

M = Dimension('mass')
L = Dimension('length')
T = Dimension('time')
I = Dimension('current')

def test_navier_stokes_case():
	""" In this case my function will be able to generate one set
	of dimensionless quantities in Navier-Stokes equations ......"""
	Lc = L # Characteristic length scale
	Uc = L/T # Characteristic velocity scale
	rho_c = M/L**3 # Characteristic density scale
	mu = M/L/T # Dynamic viscosity
	p_c = M/L/T**2 # Characteristic pressure scale
	F_bc = M*L/T**2 # Characteristic body forces scale
	list_of_quantities = {
		'Lc': Lc,
		'Uc' : Uc,
		'rho' : rho_c,
		'mu' : mu,
		'pc' : p_c,
		'Fbc' : F_bc
	}
	set_of_dimless_quant, set_of_dimless_nums =
					DimensionSystem.get_dimensionless_numbers(list_of_quantities)
	flag = True
	for dimension in set_of_dimless_nums:
		flag = flag and dimension == Dimension(1)
	try:
		assert flag
		print('\n Set of dimensionless quantities with Navier-Stokes equations '+\
			'-> \n')
		for element in set_of_dimless_quant: print(element)
	except Exception as exception: print(exception)
	""" You should be able to see the Reynolds number, Euler's
	number and Froude number in the output """
def test_mhd_case():
	""" In this case my function will be able to generate one set of
	dimensionless quantities in equations governing magnetohydrodynamics"""
	Lc = L # Characteristic length scale
	Uc = L/T # Characteristic velocity scale
	rho_c = M/L**3 # Characteristic density scale
	mu = M/L/T # Dynamic viscosity
	j0 = I/L**2 # Characteristic range value of current density
	sigma_0 = M*L**3*T**(-3)*I**(-2) # Characteristic conductivity scale
	B0 = M*T**(-2)*I**(-1) # Magnetic field characteristic scale
	list_of_quantities = {
		'Lc': Lc,
		'Uc' : Uc,
		'rho' : rho_c,
		'mu' : mu,
		'j0' : j0,
		'sigma_0' : sigma_0,
		'B0' : B0
	}
	set_of_dimless_quant, set_of_dimless_nums =
					DimensionSystem.get_dimensionless_numbers(list_of_quantities)
	flag = True
	for dimension in set_of_dimless_nums:
		flag = flag and dimension == Dimension(1)
	try:
		assert flag
		print('\n Set of dimensionless quantities with MHD equations -> \n')
		for element in set_of_dimless_quant: print(element)
	except Exception as exception: print(exception)
	""" You should be able to see the Reynolds number and the other two numbers
		contain information of the Hartmann number, Chandrasekhar number """
def test_qed_case():
	""" In this case my function will be able to generate one set
	of dimensionless quantities in circuit quantum electrodynamics ......"""
	hbar = M*L**2/T # Dimensions of the reduced Planck's constant
	omega = 1/T # Characteristic angular resonant frequency
	epsilon = I**2*T**4*M**(-1)*L**(-3) # Vaccuum permittivity
	Vm = L**3 # Volume of the resonator
	Qe = I*T # Dimensions of the electronic charge
	Uc = L/T # Velocity dimensions
	list_of_quantities = {
		'hbar': hbar,
		'Uc' : Uc,
		'omega' : omega,
		'epsilon' : epsilon,
		'Vm' : Vm,
		'Qe' : Qe
	}
	set_of_dimless_quant, set_of_dimless_nums =
					DimensionSystem.get_dimensionless_numbers(list_of_quantities)
	flag = True
	for dimension in set_of_dimless_nums:
		flag = flag and dimension == Dimension(1)
	try:
		assert flag
		print('\n Set of dimensionless quantities in circuit QED -> \n')
		for element in set_of_dimless_quant: print(element)
	except Exception as exception: print(exception)
	""" One of the values in the output shows the
	 fine-structure interaction """
