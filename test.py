import time
from sympy import *
start_time = time.time()
x,y,z,t=symbols('x y z t')
eq = (Eq(Derivative(x(t),t),x(t)*y(t)*sin(t)), Eq(Derivative(y(t),t),y(t)**2*sin(t)))
dsolve(eq)
calc=time.time()-start_time
print calc
