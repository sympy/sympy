from sympy import Ylm, Zlm, Symbol, sqrt, pi, sin, cos, cot, exp, I, S, diff, conjugate
from sympy.functions.special.spherical_harmonics import Pl, Plm, Plmcos, Ynm, Znm, Ynm_c
from sympy.utilities.pytest import XFAIL


def test_Pl():
    x = Symbol("x")
    assert Pl(0, x) == 1
    assert Pl(1, x) == x
    assert Pl(2, x) == ((3*x**2 - 1)/2).expand()
    assert Pl(3, x) == ((5*x**3 - 3*x)/2).expand()
    assert Pl(4, x) == ((35*x**4 - 30*x**2 + 3)/8).expand()
    assert Pl(5, x) == ((63*x**5 - 70*x**3 + 15*x)/8).expand()
    assert Pl(6, x) == ((231*x**6 - 315*x**4 + 105*x**2 - 5)/16).expand()


def test_Plm():
    #http://en.wikipedia.org/wiki/Legendre_function
    x = Symbol("x")
    assert Plm(0, 0, x) == 1
    assert Plm(1, -1, x) == -Plm(1, 1, x)/2
    assert Plm(1, 0, x) == x
    assert Plm(1, 1, x) == -sqrt(1 - x**2)
    assert Plm(2, -2, x) == Plm(2, 2, x)/24
    assert Plm(2, -1, x) == -Plm(2, 1, x)/6
    assert Plm(2, 0, x) == (3*x**2 - 1)/2
    assert Plm(2, 1, x) == -3*x*sqrt(1 - x**2)
    assert Plm(2, 2, x) == 3*(1 - x**2)
    assert Plm(3, -3, x) == -Plm(3, 3, x)/720
    assert Plm(3, -2, x) == Plm(3, 2, x)/120
    assert Plm(3, -1, x) == -Plm(3, 1, x)/12
    assert Plm(3, 0, x) == (5*x**3 - 3*x)/2
    assert Plm(3, 1, x).expand() == (3*(1 - 5*x**2)/2*sqrt(1 - x**2)).expand()
    assert Plm(3, 2, x) == 15*x*(1 - x**2)
    assert Plm(3, 3, x) == -15*sqrt(1 - x**2)**3


def test_Plmcos():
    #http://en.wikipedia.org/wiki/Legendre_function
    th = Symbol("th", real=True)
    assert Plmcos(0, 0, th) == 1
    assert Plmcos(1, -1, th) == sin(th)/2
    assert Plmcos(1, 0, th) == cos(th)
    assert Plmcos(1, 1, th) == -sin(th)
    assert Plmcos(2, 0, th) == (3*cos(th)**2 - 1)/2
    assert Plmcos(2, 1, th) == -3*cos(th)*sin(th)
    assert Plmcos(2, 2, th) in [3*sin(th)**2, 3*(1 - cos(th)**2)]
    assert Plmcos(3, 0, th) == (5*cos(th)**3 - 3*cos(th))/2
    assert Plmcos(3, 1, th) == -sin(th)*(15*cos(th)**2/2 - S(3)/2)
    assert Plmcos(3, 2, th) == 15*cos(th)*sin(th)**2
    assert Plmcos(3, 3, th) == -15*sin(th)**3


def test_Ylm():
    #http://en.wikipedia.org/wiki/Spherical_harmonics
    th, ph = Symbol("theta", real=True), Symbol("phi", real=True)
    assert Ylm(0, 0, th, ph) == 1/(2*sqrt(pi))
    assert Ylm(1, -1, th, ph) == S.Half*sqrt(3/(2*pi))*sin(th)*exp(-I*ph)
    assert Ylm(1, 0, th, ph) == S.Half*sqrt(3/pi)*cos(th)
    assert Ylm(1, 1, th, ph) == -S.Half*sqrt(3/(2*pi))*sin(th)*exp(I*ph)
    assert Ylm(2, 0, th, ph).expand() == (
        S(1)/4*sqrt(5/pi)*(3*cos(th)**2 - 1)).expand()
    assert Ylm(2, 1, th, ph).expand() == (
        -S.Half*sqrt(3)*sqrt(5/(2*pi))*(sin(th)*cos(th))*exp(I*ph)).expand()

    # These last 2 return the correct answer, but the answer can be simplified
    assert Ylm(2, -2, th, ph).expand() == (-sqrt(30)*exp(-2*I*ph)*
        cos(th)**2/(8*sqrt(pi)) + sqrt(30)*exp(-2*I*ph)/(8*sqrt(pi)))
    assert Ylm(2, 2, th, ph).expand() == (-sqrt(30)*exp(2*I*ph)*
        cos(th)**2/(8*sqrt(pi)) + sqrt(30)*exp(2*I*ph)/(8*sqrt(pi)))


def test_Ynm():
    th, ph = Symbol("theta", real=True), Symbol("phi", real=True)
    from sympy.abc import n,m

    assert Ynm(0, 0, th, ph).expand(func=True) == 1/(2*sqrt(pi))
    assert Ynm(1, -1, th, ph) == -exp(-2*I*ph)*Ynm(1, 1, th, ph)
    assert Ynm(1, -1, th, ph).expand(func=True) == sqrt(6)*sqrt(-cos(th)**2 + 1)*exp(-I*ph)/(4*sqrt(pi))
    assert Ynm(1, -1, th, ph).expand(func=True) == sqrt(6)*sqrt(-cos(th)**2 + 1)*exp(-I*ph)/(4*sqrt(pi))
    assert Ynm(1, 0, th, ph).expand(func=True) == sqrt(3)*cos(th)/(2*sqrt(pi))
    assert Ynm(1, 1, th, ph).expand(func=True) == -sqrt(6)*sqrt(-cos(th)**2 + 1)*exp(I*ph)/(4*sqrt(pi))
    assert Ynm(2, 0, th, ph).expand(func=True) == 3*sqrt(5)*cos(th)**2/(4*sqrt(pi)) - sqrt(5)/(4*sqrt(pi))
    assert Ynm(2, 1, th, ph).expand(func=True) == -sqrt(30)*sqrt(-cos(th)**2 + 1)*exp(I*ph)*cos(th)/(4*sqrt(pi))
    assert Ynm(2, -2, th, ph).expand(func=True) == (-sqrt(30)*exp(-2*I*ph)*cos(th)**2/(8*sqrt(pi))
                                                    + sqrt(30)*exp(-2*I*ph)/(8*sqrt(pi)))
    assert Ynm(2, 2, th, ph).expand(func=True) == (-sqrt(30)*exp(2*I*ph)*cos(th)**2/(8*sqrt(pi))
                                                   + sqrt(30)*exp(2*I*ph)/(8*sqrt(pi)))

    assert diff(Ynm(n, m, th, ph), th) == (m*cot(th)*Ynm(n, m, th, ph)
                                           + sqrt((-m + n)*(m + n + 1))*exp(-I*ph)*Ynm(n, m + 1, th, ph))
    assert diff(Ynm(n, m, th, ph), ph) == I*m*Ynm(n, m, th, ph)

    assert conjugate(Ynm(n, m, th, ph)) == (-1)**(2*m)*exp(-2*I*m*ph)*Ynm(n, m, th, ph)

    assert Ynm(n, m, -th, ph) == Ynm(n, m, th, ph)
    assert Ynm(n, m, th, -ph) == exp(-2*I*m*ph)*Ynm(n, m, th, ph)
    assert Ynm(n, -m, th, ph) == (-1)**m*exp(-2*I*m*ph)*Ynm(n, m, th, ph)


def test_Ynm_c():
    th, ph = Symbol("theta", real=True), Symbol("phi", real=True)
    from sympy.abc import n,m

    assert Ynm_c(n, m, th, ph) == (-1)**(2*m)*exp(-2*I*m*ph)*Ynm(n, m, th, ph)


def test_Zlm():
    #http://en.wikipedia.org/wiki/Solid_harmonics#List_of_lowest_functions
    th, ph = Symbol("theta", real=True), Symbol("phi", real=True)
    assert Zlm(0, 0, th, ph) == sqrt(1/(4*pi))
    assert Zlm(1, -1, th, ph) == sqrt(3/(4*pi))*sin(th)*sin(ph)
    assert Zlm(1, 0, th, ph) == sqrt(3/(4*pi))*cos(th)
    assert Zlm(1, 1, th, ph) == sqrt(3/(4*pi))*sin(th)*cos(ph)

    assert Zlm(2, -1, th, ph) == sqrt(15)*(
        cos(ph - 2*th) - cos(ph + 2*th))/(8*sqrt(pi))
    assert Zlm(2, 0, th, ph).expand() == (S(1)/4*sqrt(5/pi) *
        (3*cos(th)**2 - 1)).expand()
    assert Zlm(2, 1, th, ph) == sqrt(15)*(
        -sin(ph - 2*th) + sin(ph + 2*th))/(8*sqrt(pi))


def test_Znm():
    th, ph = Symbol("theta", real=True), Symbol("phi", real=True)
    from sympy.abc import n,m

    assert Znm(0, 0, th, ph) == Ynm(0, 0, th, ph)
    assert Znm(1, -1, th, ph) == (-sqrt(2)*I*(Ynm(1, 1, th, ph)
                                  - exp(-2*I*ph)*Ynm(1, 1, th, ph))/2)
    assert Znm(1, 0, th, ph) == Ynm(1, 0, th, ph)
    assert Znm(1, 1, th, ph) == (sqrt(2)*(Ynm(1, 1, th, ph)
                                 + exp(-2*I*ph)*Ynm(1, 1, th, ph))/2)
    assert Znm(0, 0, th, ph).expand(func=True) == 1/(2*sqrt(pi))
    assert Znm(1, -1, th, ph).expand(func=True) == (sqrt(3)*I*sqrt(-cos(th)**2 + 1)*exp(I*ph)/(4*sqrt(pi))
                                                    - sqrt(3)*I*sqrt(-cos(th)**2 + 1)*exp(-I*ph)/(4*sqrt(pi)))
    assert Znm(1, 0, th, ph).expand(func=True) == sqrt(3)*cos(th)/(2*sqrt(pi))
    assert Znm(1, 1, th, ph).expand(func=True) == (-sqrt(3)*sqrt(-cos(th)**2 + 1)*exp(I*ph)/(4*sqrt(pi))
                                                   - sqrt(3)*sqrt(-cos(th)**2 + 1)*exp(-I*ph)/(4*sqrt(pi)))
    assert Znm(2, -1, th, ph).expand(func=True) == (sqrt(15)*I*sqrt(-cos(th)**2 + 1)*exp(I*ph)*cos(th)/(4*sqrt(pi))
                                                    - sqrt(15)*I*sqrt(-cos(th)**2 + 1)*exp(-I*ph)*cos(th)/(4*sqrt(pi)))
    assert Znm(2, 0, th, ph).expand(func=True) == 3*sqrt(5)*cos(th)**2/(4*sqrt(pi)) - sqrt(5)/(4*sqrt(pi))
    assert Znm(2, 1, th, ph).expand(func=True) == (-sqrt(15)*sqrt(-cos(th)**2 + 1)*exp(I*ph)*cos(th)/(4*sqrt(pi))
                                                   - sqrt(15)*sqrt(-cos(th)**2 + 1)*exp(-I*ph)*cos(th)/(4*sqrt(pi)))
