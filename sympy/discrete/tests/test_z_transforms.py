"""
test file for z_transforms and inverse_z_transforms
"""
from sympy.discrete.z_transforms import (z_transform,
                                         inverse_z_transform,
                                         z_pairs_properties,
                                         z_pairs_prop_inverse,
                                         round_float)
from sympy.abc import z, n
from sympy.core.symbol import symbols
from sympy.core import S
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.functions.elementary.complexes import Abs

def test_z_transform():
    """
    test from z_transform table and excercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    """
    a, b, c, = symbols('a, b, c', positive=True)
    ZT = z_transform

    # included in _z_pairs_table
    assert ZT(DiracDelta(n),n,z) == 1
    assert ZT(Heaviside(n),n,z) == z/(z - 1)
    assert ZT(cos(a*n),n,z) == \
           z*(z - cos(a))/(z**2 - 2*z*cos(a) + 1)
    assert ZT(sin(a*n),n,z) == \
           z*sin(a)/(z**2 - 2*z*cos(a) + 1)
    assert ZT(cos(b*n+c),n,z) == \
           z*(z*cos(c) - cos(b - c))/(z**2 - 2*z*cos(b) + 1)
    # shift property
    assert ZT(DiracDelta(n-a),n,z) == z**(-a)
    assert ZT(DiracDelta(n+a),n,z) == z**a
    assert ZT(b*DiracDelta(n+2),n,z) == b*z**2
    assert ZT(-b*DiracDelta(n+2),n,z) == -b*z**2
    assert ZT(Heaviside(n-a),n,z) == z**-a*z/(z - 1)
    assert ZT(Heaviside(n+a),n,z) == z**a*z/(z - 1)
    assert ZT(b*Heaviside(n+a),n,z) == b*z**a*z/(z - 1)
    # multiply by n property
    assert ZT(n*Heaviside(n),n,z) == z/(z - 1)**2
    assert ZT(n**2*Heaviside(n),n,z) == z*(z + 1)/(z - 1)**3
    assert ZT(n**3*Heaviside(n),n,z) == \
           z*(z**2 + 4*z + 1)/(z - 1)**4
    # multiply by a**n property
    assert ZT(a**n*Heaviside(n),n,z) == z/(z - a)
    assert ZT(0.5**n*Heaviside(n),n,z) == z/(z - 0.5)
    assert ZT((S.One/2)**n*Heaviside(n),n,z) == 2*z/(2*z - 1)

    # time reversal property
    assert ZT(Heaviside(-n),n,z) == -1/(z - 1)
    assert ZT(-(0.5**n)*Heaviside(-n-1),n,z) == z/(z - 0.5)

    # combined properties
    assert ZT(-n*(0.5**n)*Heaviside(-n-1),n,z) == 0.5*z/(z - 0.5)**2
    assert ZT(n**2*(0.5**n)*Heaviside(n),n,z) == \
           0.5*z*(z + 0.5)/(z - 0.5)**3
    assert ZT((0.5**n)*cos(2*n)*Heaviside(n),n,z) == \
           z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25)
    assert ZT((0.5**n)*sin(2*n)*Heaviside(n),n,z) == \
           0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25)

    # cos with phase
    Fz = ZT(cos(2*n+0.25)*Heaviside(n),n,z)
    assert round_float(Fz) == \
           z*(0.968912421710645*z + 0.178246055649492)/(z**2 - 2*z*cos(2) + 1)

    return()

def test_inverse_z_transform():
    """
    test from z_transform table and excercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    """
    a, b = symbols('a, b', positive=True)
    ZT_inv = inverse_z_transform

    # included in _z_pairs_table
    assert ZT_inv(1,z,n) == DiracDelta(n)
    assert ZT_inv(z/(z - 1),z,n) == Heaviside(n)
    assert ZT_inv(z*(z - cos(2))/(z**2 - 2*z*cos(2) + 1),z,n) == \
           cos(2*n)*Heaviside(n)
    assert ZT_inv(z*sin(2)/(z**2 - 2*z*cos(2) + 1),z,n) == \
           sin(2*n)*Heaviside(n)
    assert round_float(ZT_inv(z*(z*cos(0.5) - cos(2 - 0.5))/(z**2 - 2*z*cos(2) + 1),z,n)) == \
           cos(2*n + 0.5)*Heaviside(n)

    # shift property
    assert ZT_inv(z**(-a),z,n) == DiracDelta(n - a)
    assert ZT_inv(z**a,z,n) == DiracDelta(n + a)
    assert ZT_inv(b*z**2,z,n) == b*DiracDelta(n + 2)
    assert ZT_inv(-b*z**2,z,n) == -b*DiracDelta(n + 2)

    assert ZT_inv(z**-a*z/(z - 1),z,n) == Heaviside(n - a)
    assert ZT_inv(z**a*z/(z - 1),z,n) == Heaviside(n + a)
    assert ZT_inv(b*z**a*z/(z - 1),z,n) == b*Heaviside(n+a)

    # multiply by n property nf[n] <--> -z*diff(F[z])
    assert ZT_inv(z/(z - 1)**2,z,n) == n*Heaviside(n)
    assert ZT_inv(z*(z + 1)/(z - 1)**3,z,n) == n**2*Heaviside(n)
    assert ZT_inv(z*(z**2 + 4*z + 1)/(z - 1)**4,z,n) == n**3*Heaviside(n)

    # multiply by a**n property
    assert ZT_inv(z/(z - a),z,n) == a**n*Heaviside(n)
    assert ZT_inv(z/(z - 0.5),z,n) == 0.5**n*Heaviside(n)
    assert ZT_inv(2*z/(2*z - 1),z,n) == 0.5**n*Heaviside(n)

    # time reversal property
    assert ZT_inv(-1/(z - 1),z,n) == Heaviside(-n)

    # combined properties
    assert ZT_inv(0.5*z/(z - 0.5)**2,z,n) == 0.5**n*n*Heaviside(n)
    assert ZT_inv(0.5*z*(z + 0.5)/(z - 0.5)**3,z,n) == \
           0.5**n*n**2*Heaviside(n)
    assert ZT_inv(z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),z,n) == \
           0.5**n*cos(2*n)*Heaviside(n)
    assert ZT_inv(0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25),z,n) == \
           0.5**n*sin(2*n)*Heaviside(n)

    return()

def test_z_pairs_properties():
    """
    test from z_transform properties table and excercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla-de-propiedades/
    """
    a, b, c, = symbols('a, b, c', positive=True)
    ZT_pp = z_pairs_properties

    # included in _z_pairs_table
    assert ZT_pp(DiracDelta(n),n,z) == (1, 0, True)
    assert ZT_pp(Heaviside(n),n,z) == (z/(z - 1), Abs(z) > 1, True)
    assert ZT_pp(cos(a*n),n,z) == \
           (z*(z - cos(a))/(z**2 - 2*z*cos(a) + 1), Abs(z) > 1, True)
    assert ZT_pp(sin(a*n),n,z) == \
           (z*sin(a)/(z**2 - 2*z*cos(a) + 1), Abs(z) > 1, True)
    assert ZT_pp(cos(b*n+c),n,z) == \
           (z*(z*cos(c) - cos(b - c))/(z**2 - 2*z*cos(b) + 1),

            Abs(z) > 1, True)
    # shift property
    assert ZT_pp(DiracDelta(n-a),n,z) == (z**(-a), 0, True)
    assert ZT_pp(DiracDelta(n+a),n,z) == (z**a, 0, True)
    assert ZT_pp(b*DiracDelta(n+2),n,z) == (b*z**2, 0, True)
    assert ZT_pp(-b*DiracDelta(n+2),n,z) == (-b*z**2, 0, True)
    assert ZT_pp(Heaviside(n-a),n,z) == (z**-a*z/(z - 1), Abs(z) > 1, True)
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
    assert ZT_pp(a**n*Heaviside(n),n,z) == (z/(z - a), a*Abs(z) > 1, True)
    assert ZT_pp(0.5**n*Heaviside(n),n,z) == \
           (z/(z - 0.5), 0.5*Abs(z) > 1, True)
    assert ZT_pp((S.One/2)**n*Heaviside(n),n,z) == \
           (2*z/(2*z - 1), Abs(z)/2 > 1, True)

    # time reversal property
    assert ZT_pp(Heaviside(-n),n,z) == (-1/(z - 1), Abs(z) < 1, True)
    assert ZT_pp(-(0.5**n)*Heaviside(-n-1),n,z) == \
           (z/(z - 0.5), 0.5*Abs(z) < 1, True)

    # combined properties
    assert ZT_pp(-n*(0.5**n)*Heaviside(-n-1),n,z) == \
           (0.5*z/(z - 0.5)**2, 0.5*Abs(z) < 1, True)
    assert ZT_pp(n**2*(0.5**n)*Heaviside(n),n,z) == \
           (0.5*z*(z + 0.5)/(z - 0.5)**3, 0.5*Abs(z) > 1, True)
    assert ZT_pp((0.5**n)*cos(2*n)*Heaviside(n),n,z) == \
           (z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),
            0.5*Abs(z) > 1, True)
    assert ZT_pp((0.5**n)*sin(2*n)*Heaviside(n),n,z) == \
           (0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25), 0.5*Abs(z) > 1, True)

    # cos with fase
    Fz = ZT_pp(cos(2*n+0.25)*Heaviside(n),n,z)
    assert (round_float(Fz[0]),Fz[1],Fz[2]) == \
           (z*(0.968912421710645*z + 0.178246055649492)/(z**2 - 2*z*cos(2) + 1),
            (Abs(z) > 1), True)

    return()

def test_z_pairs_prop_inverse():
    """
    test from z_transform properties table and excercises
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla-de-propiedades/
    """
    a, b = symbols('a, b', positive=True)
    ZT_pp_inv = z_pairs_prop_inverse

    # included in _z_pairs_table
    assert ZT_pp_inv(1,z,n) == (DiracDelta(n), 0, True)
    assert ZT_pp_inv(z/(z - 1),z,n) == (Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*(z - cos(2))/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (cos(2*n)*Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*sin(2)/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (sin(2*n)*Heaviside(n), Abs(z) > 1, True)
    assert ZT_pp_inv(z*(z*cos(0.5) - cos(2 - 0.5))/(z**2 - 2*z*cos(2) + 1),z,n) == \
           (cos(2*n + 0.5)*Heaviside(n), 0, 0)

    # shift property
    assert ZT_pp_inv(z**(-a),z,n) == (DiracDelta(n - a), 0, True)
    assert ZT_pp_inv(z**a,z,n) == (DiracDelta(n + a), 0, True)
    assert ZT_pp_inv(b*z**2,z,n) == (b*DiracDelta(n + 2), 0, True)
    assert ZT_pp_inv(-b*z**2,z,n) == (-b*DiracDelta(n + 2), 0, True)

    assert ZT_pp_inv(z**-a*z/(z - 1),z,n) == \
           (Heaviside(n - a), Abs(z) > 1, True)
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
           (a**n*Heaviside(n), (Abs(z)/a > 1), True)
    assert ZT_pp_inv(z/(z - 0.5),z,n) == \
           (0.5**n*Heaviside(n), (2*Abs(z) > 1), True)
    assert ZT_pp_inv(2*z/(2*z - 1),z,n) == \
           (0.5**n*Heaviside(n), 2*Abs(z) > 1, True)

    # time reversal property
    assert ZT_pp_inv(-1/(z - 1),z,n) == (Heaviside(-n), Abs(z) < 1, True)

    # combined properties
    assert ZT_pp_inv(0.5*z/(z - 0.5)**2,z,n) == \
           (0.5**n*n*Heaviside(n), 2*Abs(z) > 1, True)
    assert ZT_pp_inv(0.5*z*(z + 0.5)/(z - 0.5)**3,z,n) == \
           (0.5**n*n**2*Heaviside(n), 2*Abs(z) > 1, True)
    assert ZT_pp_inv(z*(z - 0.5*cos(2))/(z**2 - z*cos(2) + 0.25),z,n) == \
           (0.5**n*cos(2*n)*Heaviside(n), 2*Abs(z) > 1, True)
    assert ZT_pp_inv(0.5*z*sin(2)/(z**2 - z*cos(2) + 0.25),z,n) == \
           (0.5**n*sin(2*n)*Heaviside(n), 2*Abs(z) > 1, True)

    return()
