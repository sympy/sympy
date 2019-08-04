from sympy import *
x=symbols('x')
f1,f2,f3=fps(exp(x)), fps(sin(x)), fps(cos(x))
fcomp = f1.compose(f2)
print(fcomp.ak)
print(fcomp._eval_terms(8))

