import sympy as sym
x, g, K, m, t = sym.symbols('x, g, K, m, t', positive=True, real=True)
dxx = sym.Eq(x(t).diff(t, 2), g - K/m*x(t).diff(t))
ics = {x(0): 0, x(t).diff(t).subs(t, 0): 0}
print(sym.dsolve(dxx, ics=ics))
