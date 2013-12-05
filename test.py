from sympy import Poly, roots, Symbol
from sympy.utilities.solution import reset_solution

x = Symbol('x')
reset_solution()
roots(Poly((x - 2) * (x**2 + 3) * (x - 6), x))
