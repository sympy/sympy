from sympy import *
from sympy.calculus.util import is_positive_over

x, y = symbols('x y')

ex = 2**y - y
print(ex.is_positive)
print(is_positive_over(ex, y))