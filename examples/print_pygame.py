import sys
sys.path.append("..")
sys.path.append(".")

from sympy import sqrt,exp,log,Symbol,oo,Rational,sin,cos,limit,I,pi,Mul
from sympy.core import basic
from sympy.printing import print_pygame, pprint

x=Symbol("x") 
y=Symbol("y") 

expressions = (
    x**x,
    x+y+x,
    sin(x)**x,
    sin(x)**cos(x),
    sin(x)/(cos(x)**2 * x**x +(2*y)),

    sin(x**2+exp(x)),
    sqrt(exp(x)),

    #print (1/cos(x)).series(x,10)
)

print "sympy print:"
for expr in expressions:
    print expr
print ''

print "pygame print:"
for expr in expressions:
    pprint(expr)
    print_pygame(expr)
