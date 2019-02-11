from sympy import *
x =Symbol("x")
f = Function("g")(x)
g = Function("g")(x)
eq = x *f.diff(x)**2 + f.diff(x,2)
#eq = dsolve(dsolve(eq.subs(f.diff(x), g)).subs(g, f.diff(x)))

print(dsolve(eq))

