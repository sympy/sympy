from sympy import *
x=symbols('x')
f1,f2=fps(exp(x)), fps(sin(x))
fprod = f1.compose(f2)
print(fprod._eval_terms(6))

