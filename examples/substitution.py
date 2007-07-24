import sys
sys.path.append("..")

import sympy
x=sympy.Symbol('x')
y=sympy.Symbol('y')
e=1/sympy.cos(x)
print e
print e.subs(sympy.cos(x),y)
print e.subs(sympy.cos(x),y).subs(y,x**2)
e=1/sympy.log(x)
e=e.subs(x,sympy.Real("2.71828"))
print e
print e.evalf()
