from sympy.solvers.diophantine import diop_solve
from sympy import symbols

def test_diop_solve():
    x, y, t = symbols("x, y, t", type = int)

    assert diop_solve(2*x + 3*y - 5) == {x: 15*t - 5, y: -10*t + 5}
