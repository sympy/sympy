# evalf issue
from sympy import *
x = betainc(1, 2, 3, 3)
x.evalf()

# lambdify issue
# from sympy import *
# from sympy.abc import x,y,a,b
# f=lambdify([a, b], beta(a, b))
# f=lambdify([a, b, x, y], betainc(a, b, x, y))
# f(1,2)