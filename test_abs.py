from sympy import reset_solution, last_solution
from sympy import symbols
from sympy import solve
from sympy import Poly
from sympy import Abs

x, y = symbols('x, y')

print solve([Abs(x) - 3, x - 3], [x], dict = True)