from sympy.utilities.solution import reset_solution, last_solution
from sympy import symbols
from sympy.solvers.solvers import solve

reset_solution()
x, y = symbols('x,y')
solve([x + 5*y - 2, -3*x + 6*y - 15], [x, y])
R = last_solution()
for i in R:
	print i