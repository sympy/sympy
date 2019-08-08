from sympy import *
x=symbols('x')
f1,f2,f3=fps(sin(x)), fps(exp(x)), fps(cos(x))
fprod = f1.product(f2, x)
print(f2.inverse(x).truncate(8))