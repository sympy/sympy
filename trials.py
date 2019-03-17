from sympy import *
x =Symbol("x")
f = Function("f")(x)
g = Function("g")(x)
eq = x *f.diff(x)**2 + f.diff(x,2)
#print(dsolve(eq))

