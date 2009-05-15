from sympy import jn, yn, symbols, sin, cos, pi, S

def test_jn():
    z = symbols("z")
    assert jn(0, z) == sin(z)/z
    assert jn(1, z) == sin(z)/z**2 - cos(z)/z
    assert jn(2, z) == (3/z**3-1/z)*sin(z) - (3/z**2) * cos(z)
    assert jn(3, z) == (15/z**4 - 6/z**2)*sin(z) + (1/z - 15/z**3)*cos(z)
    assert jn(4, z) == (1/z + 105/z**5 - 45/z**3)*sin(z) + \
                (-105/z**4 + 10/z**2)*cos(z)
    assert jn(5, z) == (945/z**6 - 420/z**4 + 15/z**2)*sin(z) + \
                (-1/z - 945/z**5 + 105/z**3)*cos(z)
    assert jn(6, z) == (-1/z + 10395/z**7 - 4725/z**5 + 210/z**3)*sin(z) + \
                (-10395/z**6 + 1260/z**4 - 21/z**2)*cos(z)

def test_yn():
    z = symbols("z")
    assert yn(0, z) == -cos(z)/z
    assert yn(1, z) == -cos(z)/z**2-sin(z)/z
    assert yn(2, z) == -((3/z**3-1/z)*cos(z)+(3/z**2)*sin(z))

def test_sympify_yn():
    assert S(15) in yn(3, pi).atoms()
    assert yn(3, pi) == 15/pi**4 - 6/pi**2
