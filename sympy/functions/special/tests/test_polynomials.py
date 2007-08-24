import sys
sys.path.append(".")
import py
from sympy import *

x = Symbol('x')

def test_legendre():
    assert legendre(0, x) == 1
    assert legendre(1, x) == x
    assert legendre(2, x) == ((3*x**2-1)/2).expand()
    assert legendre(3, x) == ((5*x**3-3*x)/2).expand()
    assert legendre(10, -1) == 1
    assert legendre(11, -1) == -1
    assert legendre(10, 1) == 1
    assert legendre(11, 1) == 1
    assert legendre(10, 0) != 0
    assert legendre(11, 0) == 0

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

