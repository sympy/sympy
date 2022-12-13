from sympy import * 
from sympy.abc import x
exp1 = tan(x).rewrite(sin) # cot 0 = 1
print(exp1.subs(x,0))