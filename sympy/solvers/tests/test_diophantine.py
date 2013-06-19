from sympy.solvers.diophantine import diop_solve
from sympy import symbols
from sympy import Integer
x, y, z, w, t = symbols("x, y, z, w, t", integer=True)


def test_linear():

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
    assert diop_solve(x + 3*y - 4*z + w - 6) == \
           {w: 6*t, x: -6*t - 3*y + 4*z + 6, y: y, z: z}


def test_quadratic():

    #Simple Hyperbolic case: A = C = 0 and B != 0
    assert diop_solve(3*x*y + 34*x - 12*y + 1) == \
        set([(-Integer(133), -Integer(11)), (Integer(5), -Integer(57))])
    assert diop_solve(6*x*y + 2*x + 3*y + 1) == set([])
    assert diop_solve(-13*x*y + 2*x - 4*y - 54) == set([(Integer(27), Integer(0))])
    assert diop_solve(-27*x*y - 30*x - 12*y - 54) == set([(-Integer(14), -Integer(1))])
    assert diop_solve(2*x*y + 5*x + 56*y + 7) == set([(-Integer(161), -Integer(3)),\
        (-Integer(47),-Integer(6)), (-Integer(35), -Integer(12)), (-Integer(29), -Integer(69)),\
        (-Integer(27), Integer(64)), (-Integer(21), Integer(7)),(-Integer(9), Integer(1)),\
        (Integer(105), -Integer(2))])
    assert diop_solve(6*x*y + 9*x + 2*y + 3) == set([])
    assert diop_solve(x*y + x + y + 1) == set([(-Integer(1), t), (t, -Integer(1))])

    #Elliptical case: B**2 - 4AC < 0
    assert diop_solve(42*x**2 + 8*x*y + 15*y**2 + 23*x + 17*y - 4915) == set([(-Integer(11), -Integer(1))])
    assert diop_solve(4*x**2 + 3*y**2 + 5*x - 11*y + 12) == set([])
    assert diop_solve(x**2 + y**2 + 2*x + 2*y + 2) == set([(-Integer(1), -Integer(1))])
    assert diop_solve(15*x**2 - 9*x*y + 14*y**2 - 23*x - 14*y - 4950) == set([(-Integer(15), Integer(6))])
    assert diop_solve(10*x**2 + 12*x*y + 12*y**2 - 34) == set([(Integer(1), -Integer(2)),\
        (-Integer(1), -Integer(1)),(Integer(1), Integer(1)), (-Integer(1), Integer(2))])
    assert diop_solve(3*x**2 + 5*x*y + 7*y**2) == set([(Integer(0), Integer(0))])
