from sympy import cos, sin, cot, I, hyper, pprint, sec, tan, besselj, besselk, besseli, bessely
from sympy.abc import x
expr = sin(x)
print(expr.rewrite(besselj))
pprint(expr.rewrite(besselj))


