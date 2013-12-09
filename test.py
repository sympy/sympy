# from sympy import Poly, roots, Symbol
from sympy.utilities.solution import reset_solution, last_solution
from sympy import symbols
from sympy.solvers.solvers import solve
from sympy.polys import Poly

reset_solution()
x, y = symbols('x,y')
print solve([x + 5*y - 2, -3*x + 6*y - 15], [x, y])
# print solve(Poly((x - 2) * (x**2 + 3) * (x - 6), x))
R = last_solution()
for i in R:
	print i