import sys
sys.path.append("..")

import sympy
a=sympy.Symbol('a')
b=sympy.Symbol('b')
e=(a+b)**5
print e
print e.expand()
