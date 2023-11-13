from sympy import symbols, sqrt, solveset
x = symbols('x')
e = sqrt(x*(2.0 - x))*(1.0 - x)/(x*(2.0 - x))
ans = solveset(e)
print(ans)