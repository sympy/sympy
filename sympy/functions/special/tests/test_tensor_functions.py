from sympy import symbols, Dij, LeviCivita, KroneckerDelta

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

def test_kronecker_delta():
    i, j, k = symbols('i j k')
    D = KroneckerDelta
    assert D(i, i) == 1
    assert D(i, i + 1) == 0
    assert D(0, 0) == 1
    assert D(0, 1) == 0
    # assert D(i, i + k) == D(0, k)
    assert D(i + k, i + k) == 1
    assert D(i + k, i + 1 + k) == 0
    assert D(i, j).subs(dict(i=1, j=0)) == 0
    assert D(i, j).subs(dict(i=3, j=3)) == 1
