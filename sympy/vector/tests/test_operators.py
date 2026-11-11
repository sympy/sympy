from __future__ import annotations
from sympy.vector import CoordSys3D, Gradient, Divergence, Curl, VectorZero, Laplacian
from sympy.printing.repr import srepr
from sympy.vector.operators import _christoffel_symbol_2nd_kind
from sympy import zeros, Matrix, sin, cos


R = CoordSys3D('R')
s1 = R.x*R.y*R.z  # type: ignore
s2 = R.x + 3*R.y**2  # type: ignore
s3 = R.x**2 + R.y**2 + R.z**2  # type: ignore
v1 = R.x*R.i + R.z*R.z*R.j  # type: ignore
v2 = R.x*R.i + R.y*R.j + R.z*R.k  # type: ignore
v3 = R.x**2*R.i + R.y**2*R.j + R.z**2*R.k  # type: ignore


def test_Gradient():
    assert Gradient(s1) == Gradient(R.x*R.y*R.z)
    assert Gradient(s2) == Gradient(R.x + 3*R.y**2)
    assert Gradient(s1).doit() == R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    assert Gradient(s2).doit() == R.i + 6*R.y*R.j


def test_Divergence():
    assert Divergence(v1) == Divergence(R.x*R.i + R.z*R.z*R.j)
    assert Divergence(v2) == Divergence(R.x*R.i + R.y*R.j + R.z*R.k)
    assert Divergence(v1).doit() == 1
    assert Divergence(v2).doit() == 3
    # issue 22384
    Rc = CoordSys3D('R', transformation='cylindrical')
    assert Divergence(Rc.i).doit() == 1/Rc.r


def test_Curl():
    assert Curl(v1) == Curl(R.x*R.i + R.z*R.z*R.j)
    assert Curl(v2) == Curl(R.x*R.i + R.y*R.j + R.z*R.k)
    assert Curl(v1).doit() == (-2*R.z)*R.i
    assert Curl(v2).doit() == VectorZero()


def test_Laplacian():
    assert Laplacian(s3) == Laplacian(R.x**2 + R.y**2 + R.z**2)
    assert Laplacian(v3) == Laplacian(R.x**2*R.i + R.y**2*R.j + R.z**2*R.k)
    assert Laplacian(s3).doit() == 6
    assert Laplacian(v3).doit() == 2*R.i + 2*R.j + 2*R.k
    assert srepr(Laplacian(s3)) == \
            'Laplacian(Add(Pow(R.x, Integer(2)), Pow(R.y, Integer(2)), Pow(R.z, Integer(2))))'


def test_christoffel_symbol_2nd_kind_cylindrical_system():
    C = CoordSys3D("C", transformation="cylindrical")
    r, theta, z = q = C.base_scalars()
    h = C.lame_coefficients()

    Gamma_r = zeros(3)
    Gamma_theta = zeros(3)
    Gamma_z = zeros(3)
    for i in range(3):
        for j in range(3):
            Gamma_r[i, j] = _christoffel_symbol_2nd_kind(h, q, i, j, 0).doit()
            Gamma_theta[i, j] = _christoffel_symbol_2nd_kind(h, q, i, j, 1).doit()
            Gamma_z[i, j] = _christoffel_symbol_2nd_kind(h, q, i, j, 2).doit()

    assert Gamma_r == Matrix([[0, 0, 0], [0, -r, 0], [0, 0, 0]])
    assert Gamma_theta == Matrix([[0, 1/r, 0], [1/r, 0, 0], [0, 0, 0]])
    assert Gamma_z == zeros(3)


def test_christoffel_symbol_2nd_kind_spherical_system():
    S = CoordSys3D("C", transformation="spherical")
    r, theta, phi = q = S.base_scalars()
    h = S.lame_coefficients()

    Gamma_r = zeros(3)
    Gamma_theta = zeros(3)
    Gamma_phi = zeros(3)
    for i in range(3):
        for j in range(3):
            Gamma_r[i, j] = _christoffel_symbol_2nd_kind(h, q, i, j, 0).doit()
            Gamma_theta[i, j] = _christoffel_symbol_2nd_kind(h, q, i, j, 1).doit()
            Gamma_phi[i, j] = _christoffel_symbol_2nd_kind(h, q, i, j, 2).doit()

    assert Gamma_r == Matrix([[0, 0, 0], [0, -r, 0], [0, 0, -r*sin(theta)**2]])
    assert Gamma_theta == Matrix([[0, 1/r, 0], [1/r, 0, 0], [0, 0, -sin(theta)*cos(theta)]])
    assert Gamma_phi == Matrix([[0, 0, 1/r], [0, 0, cos(theta)/sin(theta)], [1/r, cos(theta)/sin(theta), 0]])
