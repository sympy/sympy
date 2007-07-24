import sys
sys.path.append("..")

import sympy
a=sympy.Symbol('a')
b=sympy.Symbol('b')
e=sympy.log((a+b)**5)
print e
e=sympy.exp(e)
print e
e=sympy.log(sympy.exp((a+b)**5))
print e
