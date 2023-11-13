from sympy import symbols, sqrt, solveset
from sympy import symbols, cos, sin, Integral, sqrt, S, Piecewise, I, Abs

# x = symbols('x')
# e = sqrt(x*(2.0 - x))*(1.0 - x)/(x*(2.0 - x))
# ans = solveset(e)
# print(ans)

x, b, z = symbols('x b z')
I1 = Integral(cos(x)/(1 - (b**2)*(sin(x))**2)**(S(3)/2), x)
I1.transform(sin(x), z)
I1.transform(sin(x), z).doit()
I1 = I1.transform(sin(x), z).doit().subs(z, sin(x))
assert I1 == Piecewise((-I*sin(x)/sqrt(b**2*sin(x)**2 - 1), Abs(b**2*sin(x)**2) > 1), (sin(x)/sqrt(-b**2*sin(x)**2 + 1), True))
