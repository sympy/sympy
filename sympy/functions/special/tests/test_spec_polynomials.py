from sympy import legendre, Symbol, hermite, chebyshevu, chebyshevt, \
        chebyshevt_root, chebyshevu_root, assoc_legendre, Rational,  \
        roots, sympify, S

x = Symbol('x')

def test_legendre():
    assert legendre(0, x) == 1
    assert legendre(1, x) == x
    assert legendre(2, x) == ((3*x**2-1)/2).expand()
    assert legendre(3, x) == ((5*x**3-3*x)/2).expand()
    assert legendre(4, x) == ((35*x**4-30*x**2+3)/8).expand()
    assert legendre(5, x) == ((63*x**5-70*x**3+15*x)/8).expand()
    assert legendre(6, x) == ((231*x**6-315*x**4+105*x**2-5)/16).expand()

    assert legendre(10, -1) == 1
    assert legendre(11, -1) == -1
    assert legendre(10, 1) == 1
    assert legendre(11, 1) == 1
    assert legendre(10, 0) != 0
    assert legendre(11, 0) == 0

    assert roots(legendre(4,x), x) == {
         (Rational(3, 7) - Rational(2, 35)*30**S.Half)**S.Half: 1,
        -(Rational(3, 7) - Rational(2, 35)*30**S.Half)**S.Half: 1,
         (Rational(3, 7) + Rational(2, 35)*30**S.Half)**S.Half: 1,
        -(Rational(3, 7) + Rational(2, 35)*30**S.Half)**S.Half: 1,
    }

def test_assoc_legendre():
    Plm=assoc_legendre
    Q=(1-x**2)**Rational(1,2)

    assert Plm(0, 0, x) ==  1
    assert Plm(1, 0, x) ==  x
    assert Plm(1, 1, x) == -Q
    assert Plm(2, 0, x) ==  (3*x**2-1)/2
    assert Plm(2, 1, x) == -3*x*Q
    assert Plm(2, 2, x) ==  3*Q**2
    assert Plm(3, 0, x) ==  (5*x**3-3*x)/2
    assert Plm(3, 1, x).expand() ==  (( 3*(1-5*x**2)/2 ).expand() * Q).expand()
    assert Plm(3, 2, x) ==  15*x * Q**2
    assert Plm(3, 3, x) == -15 * Q**3

    # negative m
    assert Plm(1,-1, x) == -Plm(1, 1, x)/2
    assert Plm(2,-2, x) ==  Plm(2, 2, x)/24
    assert Plm(2,-1, x) == -Plm(2, 1, x)/6
    assert Plm(3,-3, x) == -Plm(3, 3, x)/720
    assert Plm(3,-2, x) ==  Plm(3, 2, x)/120
    assert Plm(3,-1, x) == -Plm(3, 1, x)/12


def test_chebyshev():
    assert chebyshevt(0, x) == 1
    assert chebyshevt(1, x) == x
    assert chebyshevt(2, x) == 2*x**2-1
    assert chebyshevt(3, x) == 4*x**3-3*x
    for n in range(1, 4):
        for k in range(n):
            z = chebyshevt_root(n, k)
            assert chebyshevt(n, z) == 0
    for n in range(1, 4):
        for k in range(n):
            z = chebyshevu_root(n, k)
            assert chebyshevu(n, z) == 0

def test_hermite():
    assert hermite(6, x) == 64*x**6 - 480*x**4 + 720*x**2 - 120

