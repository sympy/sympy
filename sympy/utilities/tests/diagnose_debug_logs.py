from sympy import symbols, cos, integrate
x = symbols("x")
y = integrate((1-cos(x))/x, x)
