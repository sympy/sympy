from sympy import symbols, Dij, LeviCivita

x, y = symbols('x,y')

def test_Dij():
    assert Dij(1, 1) == 1
    assert Dij(1, 2) == 0
    assert Dij(x, x) == 1
    assert Dij(x**2-y**2, x**2-y**2) == 1

def test_levicivita():
    assert LeviCivita(1, 2, 3) == 1
    assert LeviCivita(1, 3, 2) == -1
    assert LeviCivita(1, 2, 2) == 0
    i,j,k = symbols('i j k')
    assert LeviCivita(i, j, k) == LeviCivita(i,j,k, evaluate=False)
    assert LeviCivita(i, j, i) == 0
    assert LeviCivita(1, i, i) == 0
    assert LeviCivita(i, j, k).doit() == (j - i)*(k - i)*(k - j)/2
    assert LeviCivita(1, 2, 3, 1) == 0
    assert LeviCivita(4, 5, 1, 2, 3) == 1
    assert LeviCivita(4, 5, 2, 1, 3) == -1


