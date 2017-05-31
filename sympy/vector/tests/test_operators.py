from sympy.vector import CoordSysCartesian, Gradient, Divergence, Curl, VectorZero


R = CoordSysCartesian('R')
s1 = R.x*R.y*R.z
s2 = R.x + 3*R.y**2
v1 = R.x*R.i + R.z*R.z*R.j
v2 = R.x*R.i + R.y*R.j + R.z*R.k


def test_Gradient():
    assert Gradient(s1, R) == Gradient(R.x*R.y*R.z, R)
    assert Gradient(s2, R) == Gradient(R.x + 3*R.y**2, R)
    assert Gradient(s1, R).doit() == R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    assert Gradient(s2, R).doit() == R.i + 6*R.y*R.j


def test_Divergence():
    assert Divergence(v1, R) == Divergence(R.x*R.i + R.z*R.z*R.j, R)
    assert Divergence(v2, R) == Divergence(R.x*R.i + R.y*R.j + R.z*R.k, R)
    assert Divergence(v1, R).doit() == 1
    assert Divergence(v2, R).doit() == 3


def test_Curl():
    assert Curl(v1, R) == Curl(R.x*R.i + R.z*R.z*R.j, R)
    assert Curl(v2, R) == Curl(R.x*R.i + R.y*R.j + R.z*R.k, R)
    assert Curl(v1, R).doit() == (-2*R.z)*R.i
    assert Curl(v2, R).doit() == VectorZero()
