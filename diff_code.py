from sympy import *
x=symbols('x')
e = sin(1/x + exp(-x)) - sin(1/x)
print(e.aseries(x, n=3, hir=True))

e1 = exp(exp(x)/(1 - 1/x))
print(e1.aseries(x))