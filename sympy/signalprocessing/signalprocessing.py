from sympy import fraction
from sympy.solvers import solve

def pole_zero(expression, symbol):
    num, den = fraction(expression)
    zeros = solve(num, symbol)
    poles = solve(den, symbol)

    return poles, zeros



