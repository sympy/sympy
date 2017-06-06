from sympy.physics.mechanics import MovingRefFrame, dynamicsymbols, Dyadic


A = MovingRefFrame('A')


def test_dyadic():
    d1 = A.x | A.x
    d2 = A.y | A.y
    d3 = A.x | A.y
    assert d1 * 0 == 0
    assert d1 != 0
    assert d1 * 2 == 2 * A.x | A.x
    assert d1 / 2. == 0.5 * d1
    assert d1 & (0 * d1) == 0
    assert d1 & d2 == 0
    assert d1 & A.x == A.x
    assert d1 ^ A.x == 0
    assert d1 ^ A.y == A.x | A.z
    assert d1 ^ A.z == - A.x | A.y
    assert d2 ^ A.x == - A.y | A.z
    assert A.x ^ d1 == 0
    assert A.y ^ d1 == - A.z | A.x
    assert A.z ^ d1 == A.y | A.x
    assert A.x & d1 == A.x
    assert A.y & d1 == 0
    assert A.y & d2 == A.y
    assert d1 & d3 == A.x | A.y
    assert d3 & d1 == 0
    assert d1.dt(A) == 0
    q = dynamicsymbols('q')
    qd = dynamicsymbols('q', 1)
    B = MovingRefFrame('B', orient_type = 'Axis', orient_amount = [q, A.z],
                       parentframe = A)
    assert d1.express(B) == d1.express(B, B)
    assert d1.express(B) == ((cos(q)**2) * (B.x | B.x) + (-sin(q) * cos(q)) *
            (B.x | B.y) + (-sin(q) * cos(q)) * (B.y | B.x) + (sin(q)**2) *
            (B.y | B.y))
    assert d1.express(B, A) == (cos(q)) * (B.x | A.x) + (-sin(q)) * (B.y | A.x)
    assert d1.express(A, B) == (cos(q)) * (A.x | B.x) + (-sin(q)) * (A.x | B.y)
    assert d1.dt(B) == (-qd) * (A.y | A.x) + (-qd) * (A.x | A.y)
