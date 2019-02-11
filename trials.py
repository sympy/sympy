from sympy import *
x =Symbol("x")
f = Function("f")(x)
print(dsolve(x *f.diff(x)**2 + f.diff(x,2)))
