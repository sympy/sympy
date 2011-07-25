from sympy import symbols, Dij, Eijk

x, y = symbols('x,y')

def test_Dij():
    assert Dij(1, 1) == 1
    assert Dij(1, 2) == 0
    assert Dij(x, x) == 1
    assert Dij(x**2-y**2, x**2-y**2) == 1

def test_Eijk():
    assert Eijk(1, 2, 3) == 1
    assert Eijk(2, 3, 1) == 1
    assert Eijk(3, 2, 1) == -1
    assert Eijk(1, 1, 2) == 0
    assert Eijk(1, 3, 1) == 0
    assert Eijk(1, x, x) == 0
    assert Eijk(x, y, x) == 0
