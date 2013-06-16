from sympy.solvers.diophantine import diop_solve
from sympy import symbols
from sympy import Integer


def test_linear():
    x, y, z, w, t = symbols("x, y, z, w, t", type = Integer)

    assert diop_solve(2*x + 3*y - 5) == {x: 15*t - 5, y: -10*t + 5}
    assert diop_solve(3*y + 2*x - 5) == {x: 15*t - 5, y: -10*t + 5}
    assert diop_solve(2*x - 3*y - 5) == {x: -15*t - 5, y: -10*t - 5}
    assert diop_solve(-2*x - 3*y - 5) == {x: -15*t + 5, y: 10*t - 5}
    assert diop_solve(7*x + 5*y) == {x: 5*t, y: -7*t}
    assert diop_solve(2*x + 4*y) == {x: 2*t, y: -t}
    assert diop_solve(4*x + 6*y - 4) == {x: 6*t - 2, y: -4*t + 2}
    assert diop_solve(4*x + 6*y - 3) == {x: None, y: None}
    assert diop_solve(4*x + 3*y -4*z + 5) == \
           {x: -15*t + 4*z - 5, y: 20*t - 4*z + 5, z: z}
    assert diop_solve(4*x + 2*y + 8*z - 5) == {x: None, y: None, z: None}
    assert diop_solve(5*x + 7*y - 2*z - 6) == \
           {x: 42*t + 6*z + 18, y: -30*t - 4*z - 12, z: z}
    assert diop_solve(3*x - 6*y + 12*z - 9) == \
           {x: -6*t - 4*z + 3, y: -3*t, z: z}
    assert diop_solve(x + 3*y - 4*z + w -6) == \
           {w: 6*t, x: -6*t - 3*y + 4*z + 6, y: y, z: z}
