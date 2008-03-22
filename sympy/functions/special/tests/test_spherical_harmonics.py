from sympy import Ylm, Zlm, Symbol, sympify, sqrt, pi, sin, cos, exp, I
from sympy.functions.special.spherical_harmonics import Pl, Plm, Plmcos

def test_Pl():
    x = Symbol("x")
    assert Pl(0, x) == 1
    assert Pl(1, x) == x
    assert Pl(2, x) == ((3*x**2 - 1)/2).expand()
    assert Pl(3, x) == ((5*x**3 - 3*x)/2).expand()
    assert Pl(4, x) == ((35*x**4-30*x**2+3)/8).expand()
    assert Pl(5, x) == ((63*x**5-70*x**3+15*x)/8).expand()
    assert Pl(6, x) == ((231*x**6-315*x**4+105*x**2-5)/16).expand()

def test_Plm():
    #http://en.wikipedia.org/wiki/Legendre_function
    x = Symbol("x")
    assert Plm(0, 0, x) == 1
    assert Plm(1, -1, x) == -Plm(1, 1, x)/2
    assert Plm(1, 0, x) == x
    assert Plm(1, 1, x) == -(1-x**2)**(sympify(1)/2)
    assert Plm(2, -2, x) == Plm(2, 2, x)/24
    assert Plm(2, -1, x) == -Plm(2, 1, x)/6
    assert Plm(2, 0, x) == (3*x**2-1)/2
    assert Plm(2, 1, x) == -3*x*(1-x**2)**(sympify(1)/2)
    assert Plm(2, 2, x) == 3*(1-x**2)
    assert Plm(3, -3, x) == -Plm(3, 3, x)/720
    assert Plm(3, -2, x) == Plm(3, 2, x)/120
    assert Plm(3, -1, x) == -Plm(3, 1, x)/12
    assert Plm(3, 0, x) == (5*x**3-3*x)/2
    assert Plm(3, 1, x).expand() == (( 3*(1-5*x**2)/2 ).expand() \
            *(1-x**2)**(sympify(1)/2)).expand()
    assert Plm(3, 2, x) == 15*x*(1-x**2)
    assert Plm(3, 3, x) == -15*(1-x**2)**(sympify(3)/2)

def test_Plmcos():
    #http://en.wikipedia.org/wiki/Legendre_function
    th = Symbol("th", real = True)
    assert Plmcos(0, 0, th) == 1
    assert Plmcos(1, -1, th) == sin(th)/2
    assert Plmcos(1, 0, th) == cos(th)
    assert Plmcos(1, 1, th) == -sin(th)
    assert Plmcos(2, 0, th) == (3*cos(th)**2-1)/2
    assert Plmcos(2, 1, th) == -3*cos(th)*sin(th)
    assert Plmcos(2, 2, th) in [3*sin(th)**2, 3*(1-cos(th)**2)]
    assert Plmcos(3, 0, th) == (5*cos(th)**3-3*cos(th))/2
    assert Plmcos(3, 1, th) == -3*(5*cos(th)**2-1)/2 *sin(th)
    assert Plmcos(3, 2, th) == 15*cos(th)*sin(th)**2
    assert Plmcos(3, 3, th) == -15*sin(th)**3

def test_Ylm():
    #http://en.wikipedia.org/wiki/Spherical_harmonics
    th, ph = Symbol("theta", real = True), Symbol("phi", real = True)
    assert Ylm(0, 0, th, ph) == sympify(1)/(2*sqrt(pi))
    assert Ylm(1, -1, th, ph) == sympify(1)/2 * sqrt(3/(2*pi)) * sin(th) * \
            exp(-I*ph)
    assert Ylm(1, 0, th, ph) == sympify(1)/2 * sqrt(3/pi) * cos(th)
    assert Ylm(1, 1, th, ph) == -sympify(1)/2 * sqrt(3/(2*pi)) * sin(th) * \
            exp(I*ph)
    #Ylm returns here a correct, but different expression:
    #assert Ylm(2, -2, th, ph).expand() == (sympify(1)/4 * sqrt(15/(2*pi)) * \
    #        sin(th)**2 * exp(-2*I*ph)).expand()
    assert Ylm(2, 0, th, ph).expand() == (sympify(1)/4 * sqrt(5/pi) * \
            (3*cos(th)**2-1)).expand()
    assert Ylm(2, 1, th, ph).expand() == (-sympify(1)/2 * \
            sqrt(3)*sqrt(5/(2*pi)) * (sin(th)*cos(th)) * exp(I*ph)).expand()
    #Ylm returns here a correct, but different expression:
    #assert Ylm(2, 2, th, ph).expand() == (sympify(1)/4 * sqrt(15/(2*pi)) * \
    #        sin(th)**2 * exp(2*I*ph)).expand()

def test_Zlm():
    #http://en.wikipedia.org/wiki/Solid_harmonics#List_of_lowest_functions
    th, ph = Symbol("theta", real = True), Symbol("phi", real = True)
    assert Zlm(0, 0, th, ph) == sqrt(1/(4*pi))
    assert Zlm(1, -1, th, ph) == sqrt(3/(4*pi))*sin(th)*sin(ph)
    assert Zlm(1, 0, th, ph) == sqrt(3/(4*pi))*cos(th)
    assert Zlm(1, 1, th, ph) == sqrt(3/(4*pi))*sin(th)*cos(ph)

    assert Zlm(2, -1, th, ph) == sqrt(15/(4*pi))*sin(th)*cos(th)*sin(ph)
    assert Zlm(2, 0, th, ph).expand() == (sympify(1)/4 * sqrt(5/pi) * \
            (3*cos(th)**2-1)).expand()
    assert Zlm(2, 1, th, ph) == sqrt(15/(4*pi))*sin(th)*cos(th)*cos(ph)
