import sys
sys.path.append("..")
sys.path.append(".")

from sympy import exp,log,Symbol,oo,Rational,sin,cos,limit,I,pi,Mul
from sympy.core import basic

x=Symbol("x") 
y=Symbol("y") 

def p():
    print x**x
    print x+y+x
    print sin(x)**x
    print sin(x)**cos(x)
    print sin(x)/(cos(x)**2 * x**x +(2*y))

    print sin(x**2+exp(x))
    print sqrt(exp(x))
    print sqrt(exp(x))

    #print (1/cos(x)).series(x,10)

print "sympy print:"
p()
#basic.outputType="tex"
basic.outputType="pygame"
print "_"*70
print "pretty print:"
p()
