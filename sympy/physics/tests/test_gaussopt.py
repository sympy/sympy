from sympy import N, Symbol, Matrix, symbols, factor, I, pi, Float, sqrt
from sympy.physics.gaussopt import (BeamParameter, CurvedMirror, 
  CurvedRefraction, FlatMirror, FlatRefraction, FreeSpace, GeometricRay, 
  FreeSpace, RayTransferMatrix, ThinLens, ThinLens, conjugate_gauss_beams, 
  gaussian_conj , geometric_conj_ab, geometric_conj_af, geometric_conj_bf, 
  rayleigh2waist, waist2rayleigh)
def streq(a, b):
    return str(a) == str(b)
def test_gauss_opt():
    mat = RayTransferMatrix(1,2,3,4)
    assert mat == Matrix([[1, 2],[3, 4]])
    assert mat == RayTransferMatrix( Matrix([[1,2],[3,4]]) )
    assert [mat.A, mat.B, mat.C, mat.D] == [1, 2, 3, 4]
    f = Symbol('f')
    lens = ThinLens(f)
    assert lens == Matrix([[   1, 0], [-1/f, 1]])
    assert lens.C == -1/f
    d = symbols('d')
    assert FreeSpace(d) == Matrix([[ 1, d], [0, 1]])
    n1, n2 = symbols('n1 n2')
    assert FlatRefraction(n1, n2) == Matrix([[1,     0], [0, n1/n2]])
    R, n1, n2 = symbols('R n1 n2')
    assert CurvedRefraction(R, n1, n2) == Matrix([[1, 0], [(n1 - n2)/(R*n2), n1/n2]])
    assert FlatMirror() == Matrix([[1, 0], [0, 1]])
    R = symbols('R')
    assert CurvedMirror(R) == Matrix([[   1, 0], [-2/R, 1]])
    f = symbols('f')
    assert ThinLens(f) == Matrix([[   1, 0], [-1/f, 1]])
    d,h,angle = symbols('d,h,angle')
    assert GeometricRay(h,angle) == Matrix([[    h], [angle]])
    assert FreeSpace(d)*GeometricRay(h,angle) == Matrix([[angle*d + h], [angle]])
    assert GeometricRay( Matrix( ((h,),(angle,)) ) ) == Matrix([[h], [angle]])
    p = BeamParameter(530e-9, 1, w=1e-3)
    assert streq(p.q, 1 + 1.88679245283019*I*pi)
    assert streq(N(p.q), 1.0 + 5.92753330865999*I)
    assert streq(N(p.w_0), Float(0.00100000000000000))
    assert streq(N(p.z_r), Float(5.92753330865999))
    fs = FreeSpace(10)
    p1 = fs*p
    assert streq(N(p.w), Float(0.00101413072159615))
    assert streq(N(p1.w), Float(0.00210803120913829))
    w, wavelen = symbols('w wavelen')
    assert waist2rayleigh(w, wavelen) == pi*w**2/wavelen
    z_r, wavelen = symbols('z_r wavelen')
    assert rayleigh2waist(z_r, wavelen) == sqrt(wavelen*z_r)/sqrt(pi)
    a, b = symbols('a b')
    assert geometric_conj_ab(a, b) == a*b/(a + b)
    a, b, f = symbols('a b f')
    assert geometric_conj_af(a, f) == a*f/(a - f)
    assert geometric_conj_bf(b, f) == b*f/(b - f)
    s_in, z_r_in, f = symbols('s_in z_r_in f')
    assert gaussian_conj(s_in, z_r_in, f)[0] == 1/(-1/(s_in + z_r_in**2/(-f + s_in)) + 1/f)
    assert gaussian_conj(s_in, z_r_in, f)[1] == z_r_in/(1 - s_in**2/f**2 + z_r_in**2/f**2)
    assert gaussian_conj(s_in, z_r_in, f)[2] == 1/sqrt(1 - s_in**2/f**2 + z_r_in**2/f**2)
    l, w_i, w_o, f = symbols('l w_i w_o f')
    assert conjugate_gauss_beams(l, w_i, w_o, f=f)[0] == f*(-sqrt(w_i**2/w_o**2 - pi**2*w_i**4/(f**2*l**2)) + 1)
    assert factor(conjugate_gauss_beams(l, w_i, w_o, f=f)[1]) == f*w_o**2*(w_i**2/w_o**2 - sqrt(w_i**2/w_o**2 - pi**2*w_i**4/(f**2*l**2)))/w_i**2
    assert conjugate_gauss_beams(l, w_i, w_o, f=f)[2] == f
