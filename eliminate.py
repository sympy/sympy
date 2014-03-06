from sympy.core.relational import Relational

	
def eliminate(eq1,eq2,k):
	'''
	This function eliminates 'k' variable from equation 2 and returns equation 1 without the 'k' term.
	Parameters:
	eq1 : Equation 1
	eq2 : Equation 2
	k   : Variable to be eliminated
	
	Examples
    ========
	>>>eliminate(x**2 + y**2 + z , y +z,z)
	x2+y2=y
	
	>>>eliminate(x**2 + y**2 + z , y +z + x*y,z)
	x2+y2=y(x+1)

	>>>eliminate(x**2 + y**3 +z, x**3 + z, z)
	x2+y3=x3


	'''
	eqq1 = first_eq(eq1)
	eqq2 = second_eq(eq2)
	elim_var = k
	elim_var2 = solve(eqq2,elim_var)
	elim_var1 = solve(eqq1,elim_var)
	return Relational(-elim_var1[0],-elim_var2[0],'==')
	
'''==============================================================================================================='''
	
def first_eq(eq):
	return Eq(eq)
	
def second_eq(eq):
	return Eq(eq)
	
'''==============================================================================================================='''
