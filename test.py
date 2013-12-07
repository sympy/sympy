from sympy import Poly, roots, Symbol
from sympy.utilities.solution import reset_solution, last_solution

x = Symbol('x')
reset_solution()
print roots(Poly(6 * x ** 3 - 11 * x ** 2 - 2 * x + 8, x))
R = last_solution()
for i in R:
	print i