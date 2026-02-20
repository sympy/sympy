"""
test file for z_transforms and inverse_z_transforms
2023/04/27
"""
from sympy.integrals.ztransforms import (z_transform,
                                         inverse_z_transform,
                                         _z_pairs_properties,
                                         _z_pairs_prop_inverse,
                                         _round_float)
from sympy.abc import z, n
#from sympy.core.numbers import Rational
from sympy.core.symbol import symbols
from sympy.core import S
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.functions.elementary.complexes import Abs

def test_z_transform():
    """
    test from z_transform table and exercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    """
    a, b, c, = symbols('a, b, c', positive=True)
    ZT = z_transform

    # included in _z_pairs_table
    assert ZT(DiracDelta(n),n,z) == (1, True, True)
    assert ZT(Heaviside(n),n,z) == (z/(z - 1), Abs(z) > 1, True)
    assert ZT(cos(a*n),n,z) == \
           (z*(z - cos(a))/(z**2 - 2*z*cos(a) + 1), Abs(z) > 1, True)
    assert ZT(sin(a*n),n,z) == \
           (z*sin(a)/(z**2 - 2*z*cos(a) + 1), Abs(z) > 1, True)
    assert ZT(cos(b*n+c),n,z) == \
           (z*(z*cos(c) - cos(b - c))/(z**2 - 2*z*cos(b) + 1),
            Abs(z) > 1, True)
    # shift property
    assert ZT(DiracDelta(n-a),n,z) == (z**(-a), Abs(z) > 0, True)
    assert ZT(DiracDelta(n+a),n,z) == (z**a, True, True)
    assert ZT(b*DiracDelta(n+2),n,z) == (b*z**2, True, True)
    assert ZT(-b*DiracDelta(n+2),n,z) == (-b*z**2, True, True)
    assert ZT(Heaviside(n-a),n,z) == (z**(1-a)/(z - 1), Abs(z) > 1, True)
    assert ZT(Heaviside(n+a),n,z) == (z**a*z/(z - 1), Abs(z) > 1, True)
    assert ZT(b*Heaviside(n+a),n,z) == (b*z**a*z/(z - 1), Abs(z) > 1, True)
    # multiply by n property
    assert ZT(n*Heaviside(n),n,z) == (z/(z - 1)**2, Abs(z) > 1, True)
    assert ZT(n**2*Heaviside(n),n,z) == \
            (z*(z + 1)/(z - 1)**3, Abs(z) > 1, True)
    assert ZT(n**3*Heaviside(n),n,z) == \
           (z*(z**2 + 4*z + 1)/(z - 1)**4, Abs(z) > 1, True)
    # multiply by a**n property
    assert ZT(a**n*Heaviside(n),n,z) == (z/(-a + z), a < Abs(z), True)
    assert ZT((S.One/2)**n*Heaviside(n),n,z) == \
            (z/(z - S.One/2), Abs(z) > S.One/2, True)
    assert ZT((S.Half)**n*Heaviside(n),n,z) == \
            (z/(z - S.Half), Abs(z) > S.Half, True)
    assert ZT(0.5**n*Heaviside(n),n,z) == \
            (z/(z - 0.5), Abs(z) > 0.5, True)

    # time reversal property
    assert ZT(Heaviside(-n),n,z) == (-1/(z - 1), Abs(z) < 1, True)
    assert ZT(-(S.Half**n)*Heaviside(-n-1),n,z) == \
            (z/(z - S.Half), Abs(z) < S.Half, True)
    assert ZT(-((S.One/2)**n)*Heaviside(-n-1),n,z) == \
            (z/(z - S.One/2), (Abs(z) < S.One/2), True)

    # combined properties
    assert ZT(-n*(S.Half**n)*Heaviside(-n-1),n,z) == \
            (z/(2*(z - S.Half)**2), Abs(z) < S.Half, True)
    assert ZT(-n*((S.One/2)**n)*Heaviside(-n-1),n,z) == \
              (z/(2*(z - S.One/2)**2),(Abs(z) < S.One/2), True)
    assert ZT(n**2*((S.One/2)**n)*Heaviside(n),n,z) == \
              (z*(z + S.One/2)/(2*(z - S.One/2)**3), Abs(z) > S.One/2, True)
    assert ZT(n**2*(0.5**n)*Heaviside(n),n,z) == \
              (0.5*z*(z + 0.5)/(z - 0.5)**3, Abs(z) > 0.5, True)
    assert ZT(((S.One/2)**n)*cos(2*n)*Heaviside(n),n,z) == \
              (z*(z - S.Half*cos(2))/(z**2 - z*cos(2) + S.One/4),
               Abs(z) > S.Half, True)
    assert ZT((0.5**n)*cos(2*n)*Heaviside(n),n,z) == \
              (z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),
               Abs(z) > 0.5, True)
    assert ZT(((S.One/2)**n)*sin(2*n)*Heaviside(n),n,z) == \
           (S.Half*z*sin(2)/(z**2 - z*cos(2) + S.One/4), Abs(z) > S.One/2, True)
    assert ZT((0.5**n)*sin(2*n)*Heaviside(n),n,z) == \
              (0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25), Abs(z) > 0.5, True)

    # cos with phase
    Fz = ZT(cos(2*n+0.25)*Heaviside(n),n,z)
    assert (_round_float(Fz[0]),Fz[1],Fz[2]) == \
           (z*(z*cos(0.25) - cos(1.75,evaluate=False))/(z**2 - 2*z*cos(2) + 1),
            Abs(z) > 1, True)
    return()

def test_inverse_z_transform():
    """
    test from z_transform table and exercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    """
    a, b = symbols('a, b', positive=True)
    ZT_inv = inverse_z_transform

    # included in _z_pairs_table
    assert ZT_inv(1,z,n) == (DiracDelta(n), True, True)
    assert ZT_inv(z/(z - 1),z,n) == (Heaviside(n), Abs(z) > 1, True)
    assert ZT_inv(z*(z - cos(2))/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (cos(2*n)*Heaviside(n), Abs(z) > 1, True)
    assert ZT_inv(z*sin(2)/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (sin(2*n)*Heaviside(n), Abs(z) > 1, True)
    assert ZT_inv(z*(z*cos(0.5) - cos(2 - 0.5))/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (cos(2*n + 0.5)*Heaviside(n), Abs(z) > 1, True)

    # shift property
    assert ZT_inv(z**(-a),z,n) == (DiracDelta(n-a), True, True)
    assert ZT_inv(z**a,z,n) == (DiracDelta(a + n), True, True)
    assert ZT_inv(b*z**2,z,n) == (b*DiracDelta(n + 2), True, True)
    assert ZT_inv(-b*z**2,z,n) == (-b*DiracDelta(n + 2), True, True)

    #assert ZT_inv((z**(1 - a)/(z - 1),z,n) == (Heaviside(n - a), Abs(z) > 1, True)
    assert ZT_inv(z**a*z/(z - 1),z,n) == (Heaviside(n + a), Abs(z) > 1, True)
    assert ZT_inv(b*z**a*z/(z - 1),z,n) == (b*Heaviside(n + a), Abs(z) > 1, True)

    # multiply by n property nf[n] <--> -z*diff(F[z])
    assert ZT_inv(z/(z - 1)**2,z,n) == (n*Heaviside(n), Abs(z) > 1, True)
    assert ZT_inv(z*(z + 1)/(z - 1)**3,z,n) == (n**2*Heaviside(n), Abs(z) > 1, True)
    assert ZT_inv(z*(z**2 + 4*z + 1)/(z - 1)**4,z,n) == (n**3*Heaviside(n), Abs(z) > 1, True)

    # multiply by a**n property
    assert ZT_inv(z/(z - a),z,n) == (a**n*Heaviside(n), Abs(z) > a, True)
    assert ZT_inv(z/(z - S.One/2),z,n) == (Heaviside(n)/2**n, Abs(z) > S.One/2, True)
    assert ZT_inv(z/(z - S.Half),z,n) == (Heaviside(n)/2**n, Abs(z) > S.One/2, True)
    assert ZT_inv(z/(z - 0.5),z,n) == (0.5**n*Heaviside(n), Abs(z) > 0.5, True)


    # time reversal property
    assert ZT_inv(-1/(z - 1),z,n) == (Heaviside(-n), Abs(z) < 1, True)

    # combined properties
    assert ZT_inv(0.5*z/(z - 0.5)**2,z,n) == \
           (0.5**n*n*Heaviside(n), Abs(z) > 0.5, True)
    assert ZT_inv(0.5*z*(z + 0.5)/(z - 0.5)**3,z,n) == \
           (0.5**n*n**2*Heaviside(n), Abs(z) > 0.5, True)
    assert ZT_inv(z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),z,n) == \
           (0.5**n*cos(2*n)*Heaviside(n), Abs(z) > 0.5, True)
    assert ZT_inv(0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25),z,n) == \
           (0.5**n*sin(2*n)*Heaviside(n), Abs(z) > 0.5, True)
    # Hsu, Schaum 4.24. Signal and Systems, p194
    assert ZT_inv((z/(z-3))**2,z,n) == \
           3**n*(n + 1)*Heaviside(n)

    return()

def test_z_pairs_properties():
    """
    test from z_transform properties table and exercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla-de-propiedades/
    """
    a, b, c, = symbols('a, b, c', positive=True)
    ZT_pp = _z_pairs_properties

    # included in _z_pairs_table
    assert ZT_pp(DiracDelta(n),n,z) == (1, True, True)
    assert ZT_pp(Heaviside(n),n,z) == (z/(z - 1), Abs(z) > 1, True)
    assert ZT_pp(cos(a*n),n,z) == \
           (z*(z - cos(a))/(z**2 - 2*z*cos(a) + 1), Abs(z) > 1, True)
    assert ZT_pp(sin(a*n),n,z) == \
           (z*sin(a)/(z**2 - 2*z*cos(a) + 1), Abs(z) > 1, True)
    assert ZT_pp(cos(b*n+c),n,z) == \
           (z*(z*cos(c) - cos(b - c))/(z**2 - 2*z*cos(b) + 1),
            Abs(z) > 1, True)
    # shift property
    assert ZT_pp(DiracDelta(n-a),n,z) == (z**(-a), Abs(z) > 0, True)
    assert ZT_pp(DiracDelta(n+a),n,z) == (z**a, True, True)
    assert ZT_pp(b*DiracDelta(n+2),n,z) == (b*z**2, True, True)
    assert ZT_pp(-b*DiracDelta(n+2),n,z) == (-b*z**2, True, True)
    #assert ZT_pp(Heaviside(n-a),n,z) == ((1/z**a)*z/(z - 1), Abs(z) > 1, True)
    assert ZT_pp(Heaviside(n+a),n,z) == (z**a*z/(z - 1), Abs(z) > 1, True)
    assert ZT_pp(b*Heaviside(n+a),n,z) == \
           (b*z**a*z/(z - 1), Abs(z) > 1, True)

    # multiply by n property
    assert ZT_pp(n*Heaviside(n),n,z) == (z/(z - 1)**2, Abs(z) > 1, True)
    assert ZT_pp(n**2*Heaviside(n),n,z) == \
           (z*(z + 1)/(z - 1)**3, Abs(z) > 1, True)
    assert ZT_pp(n**3*Heaviside(n),n,z) == \
           (z*(z**2 + 4*z + 1)/(z - 1)**4, Abs(z) > 1, True)

    # multiply by a**n property
    assert ZT_pp(a**n*Heaviside(n),n,z) == (z/(z - a), a < Abs(z), True)
    assert ZT_pp((S.One/2)**n*Heaviside(n),n,z) == \
           (z/(z - S.One/2), (Abs(z) > S.One/2), True)
    assert ZT_pp((S.Half)**n*Heaviside(n),n,z) == \
           (z/(z - S.Half), (Abs(z) > S.Half), True)
    assert ZT_pp(0.5**n*Heaviside(n),n,z) == \
           (z/(z - 0.5), (Abs(z) > 0.5), True)

    # time reversal property
    assert ZT_pp(Heaviside(-n),n,z) == (-1/(z - 1), Abs(z) < 1, True)
    assert ZT_pp(-(S.Half**n)*Heaviside(-n-1),n,z) == \
           (z/(z - S.Half), (Abs(z) < S.Half), True)
    assert ZT_pp(-((S.One/2)**n)*Heaviside(-n-1),n,z) == \
           (z/(z - S.One/2), (Abs(z) < S.One/2), True)
    assert ZT_pp(-((0.5)**n)*Heaviside(-n-1),n,z) == \
           (z/(z - 0.5), (Abs(z) < 0.5), True)

    # combined properties
    assert ZT_pp(-n*(S.Half**n)*Heaviside(-n-1),n,z) == \
           (z/(2*(z - S.One/2)**2),(Abs(z) < S.One/2), True)
    assert ZT_pp(-n*((S.One/2)**n)*Heaviside(-n-1),n,z) == \
           (z/(2*(z - S.One/2)**2),(Abs(z) < S.One/2), True)
    assert ZT_pp(-n*(0.5**n)*Heaviside(-n-1),n,z) == \
           (0.5*z/(z - 0.5)**2, (Abs(z) < 0.5), True)
    assert ZT_pp(n**2*((S.One/2)**n)*Heaviside(n),n,z) == \
           (z*(z + S.One/2)/(2*(z - S.One/2)**3), Abs(z) > S.One/2, True)
    assert ZT_pp(n**2*(0.5**n)*Heaviside(n),n,z) == \
           (0.5*z*(z + 0.5)/(z - 0.5)**3, Abs(z) > 0.5, True)
    assert ZT_pp(((S.One/2)**n)*cos(2*n)*Heaviside(n),n,z) == \
           (z*(z - cos(2)/2)/(z**2 - z*cos(2) + S.One/4),
            Abs(z) > S.One/2, True)
    assert ZT_pp((0.5**n)*cos(2*n)*Heaviside(n),n,z) == \
           (z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),
            Abs(z) > 0.5, True)
    assert ZT_pp(((S.One/2)**n)*sin(2*n)*Heaviside(n),n,z) == \
           (S.Half*z*sin(2)/(z**2 - z*cos(2) + S.One/4), Abs(z) > S.One/2, True)
    assert ZT_pp((0.5**n)*sin(2*n)*Heaviside(n),n,z) == \
           (0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25), Abs(z) > 0.5, True)

    # cos with fase
    Fz = ZT_pp(cos(2*n+0.25)*Heaviside(n),n,z)
    assert (_round_float(Fz[0]),Fz[1],Fz[2]) == \
            (z*(z*cos(0.25) - cos(1.75,evaluate=False))/(z**2 - 2*z*cos(2) + 1),
             Abs(z) > 1, True)
    return()

def test_z_pairs_prop_inverse():
    """
    test from z_transform properties table and exercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla-de-propiedades/
    """
    a, b = symbols('a, b', positive=True)
    ZT_pp_inv = _z_pairs_prop_inverse

    # included in _z_pairs_table
    assert ZT_pp_inv(1,z,n) == (DiracDelta(n), True, True)
    assert ZT_pp_inv(z/(z - 1),z,n) == (Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*(z - cos(2))/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (cos(2*n)*Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*sin(2)/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (sin(2*n)*Heaviside(n), Abs(z) > 1, True)
    F = ZT_pp_inv(z*(z*cos(0.5) - cos(2 - 0.5))/(z**2 - 2*z*cos(2) + 1),z,n)
    assert (_round_float(F[0]),F[1],F[2]) == \
           (cos(2*n + 0.5)*Heaviside(n), (Abs(z) > S.One), True)

    # shift property
    assert ZT_pp_inv(z**(-a),z,n) == (DiracDelta(n - a), True, True)
    assert ZT_pp_inv(z**a,z,n) == (DiracDelta(n + a), True, True)
    assert ZT_pp_inv(b*z**2,z,n) == (b*DiracDelta(n + 2), True, True)
    assert ZT_pp_inv(-b*z**2,z,n) == (-b*DiracDelta(n + 2), True, True)

    #assert ZT_pp_inv((z**(1-a)/(z - 1),z,n) == \
    #       (Heaviside(n - a), Abs(z) > 1, True)
    assert ZT_pp_inv(z**a*z/(z - 1),z,n) == \
           (Heaviside(n + a), Abs(z) > 1, True)
    assert ZT_pp_inv(b*z**a*z/(z - 1),z,n) == \
           (b*Heaviside(n+a), Abs(z) > 1, True)

    # multiply by n property
    assert ZT_pp_inv(z/(z - 1)**2,z,n) == (n*Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*(z + 1)/(z - 1)**3,z,n) == \
           (n**2*Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*(z**2 + 4*z + 1)/(z - 1)**4,z,n) == \
           (n**3*Heaviside(n), Abs(z) > 1, True)

    # multiply by a**n property
    assert ZT_pp_inv(z/(z - a),z,n) == \
           (a**n*Heaviside(n), Abs(z) > a, True)
    assert ZT_pp_inv(z/(z - S.One/2),z,n) == \
           (Heaviside(n)/2**n, Abs(z) > S.One/2, True)
    assert ZT_pp_inv(z/(z - S.Half),z,n) == \
           (Heaviside(n)/2**n, Abs(z) > S.One/2, True)
    assert ZT_pp_inv(z/(z - 0.5),z,n) == \
           (0.5**n*Heaviside(n), Abs(z) > 0.5, True)

    # time reversal property
    assert ZT_pp_inv(-1/(z - 1),z,n) == (Heaviside(-n), Abs(z) < 1, True)

    # combined properties
    assert ZT_pp_inv(0.5*z/(z - 0.5)**2,z,n) == \
           (0.5**n*n*Heaviside(n), Abs(z) > 0.5, True)
    assert ZT_pp_inv(0.5*z*(z + 0.5)/(z - 0.5)**3,z,n) == \
           (0.5**n*n**2*Heaviside(n), Abs(z) > 0.5, True)
    assert ZT_pp_inv(z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),z,n) == \
           (0.5**n*cos(2*n)*Heaviside(n), Abs(z) > 0.5, True)
    assert ZT_pp_inv(0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25),z,n) == \
           (0.5**n*sin(2*n)*Heaviside(n), Abs(z) > 0.5, True)

    return()
