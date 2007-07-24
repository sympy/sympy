import sys
sys.path.append("..")

from sympy import Symbol, cos, sin
x=Symbol('x')
e=1/cos(x)
print e.series(x,10)
e=1/sin(x)
print e.series(x,4)
