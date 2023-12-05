from sympy import symbols, sqrt, solveset
from sympy import symbols, cos, sin, Integral, sqrt, S, Piecewise, I, Abs

# right now the bug is that it is S.EmptySet when it should have 1 as a solution
x = symbols('x')
e = sqrt(x*(2.0 - x))*(1.0 - x)/(x*(2.0 - x))
ans = solveset(e)
print(ans)

# below has 1 as a solution if not float
e = sqrt(x*(2 - x))*(1 - x)/(x*(2 - x))
ans = solveset(e)
print(ans)

