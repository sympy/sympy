import sys
sys.path.append("..")

import sympy
a=sympy.Symbol('a')
b=sympy.Symbol('b')
c=sympy.Symbol('c')
e=( a*b*b+2*b*a*b )**c
print e
