import sys
sys.path.append("..")

from sympy import Symbol, cos, sin
x=Symbol('x')

e=1/cos(x)
print "Series for sec(x):"
print e.series(x,10)
print ""

e=1/sin(x)
print "Series for csc(x):"
print e.series(x,4)
