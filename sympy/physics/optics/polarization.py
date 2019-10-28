# -*- coding: utf-8 -*-
r"""
The module implements:

- Jones vectors.

- Stokes vectors.

- Jones matrices.

- Mueller matrices.

Please see the description of the individual functions for further
details and examples.

Credits and Copyright
~~~~~~~~~~~~~~~~~~~~~

AUTHOR:

- Oscar Gerardo Lazo Arjona (2017-06-18): added Wigner D matrices

Copyright (C) 2019 Oscar Gerardo Lazo Arjona <algebraicamente@gmail.com>

"""

from sympy import sin, cos, exp, I, pi, sqrt, Matrix, Abs, re, im, simplify
from sympy.physics.quantum import TensorProduct


def jones_vector(psi, chi):
    u"""A Jones vector corresponding to a polarization ellipse with `psi` tilt,
    and `chi` circularity.

    INPUT:

    -  ``psi`` - The tilt of the polarization relative to the `x` axis.
    -  ``chi`` - The angle adjacent to the mayor axis of the polarization
        ellipse.

    OUTPUT:

    A Jones vector.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> psi, chi = symbols("psi, chi", real=True)

    A general Jones vector.
    >>> pprint(jones_vector(psi, chi), use_unicode=True)
    ⎡-ⅈ⋅sin(χ)⋅sin(ψ) + cos(χ)⋅cos(ψ)⎤
    ⎢                                ⎥
    ⎣ⅈ⋅sin(χ)⋅cos(ψ) + sin(ψ)⋅cos(χ) ⎦

    Horizontal polarization
    >>> pprint(jones_vector(0, 0), use_unicode=True)
    ⎡1⎤
    ⎢ ⎥
    ⎣0⎦

    Vertical polarization
    >>> pprint(jones_vector(pi/2, 0), use_unicode=True)
    ⎡0⎤
    ⎢ ⎥
    ⎣1⎦

    Diagonal polarization
    >>> pprint(jones_vector(pi/4, 0), use_unicode=True)
    ⎡√2⎤
    ⎢──⎥
    ⎢2 ⎥
    ⎢  ⎥
    ⎢√2⎥
    ⎢──⎥
    ⎣2 ⎦

    Anti-diagonal polarization
    >>> pprint(jones_vector(-pi/4, 0), use_unicode=True)
    ⎡ √2 ⎤
    ⎢ ── ⎥
    ⎢ 2  ⎥
    ⎢    ⎥
    ⎢-√2 ⎥
    ⎢────⎥
    ⎣ 2  ⎦

    Right-hand circular polarization
    >>> pprint(jones_vector(0, pi/4), use_unicode=True)
    ⎡ √2 ⎤
    ⎢ ── ⎥
    ⎢ 2  ⎥
    ⎢    ⎥
    ⎢√2⋅ⅈ⎥
    ⎢────⎥
    ⎣ 2  ⎦

    Left-hand circular polarization
    >>> pprint(jones_vector(0, -pi/4), use_unicode=True)
    ⎡  √2  ⎤
    ⎢  ──  ⎥
    ⎢  2   ⎥
    ⎢      ⎥
    ⎢-√2⋅ⅈ ⎥
    ⎢──────⎥
    ⎣  2   ⎦

    """
    return Matrix([-I*sin(chi)*sin(psi) + cos(chi)*cos(psi),
                   I*sin(chi)*cos(psi) + sin(psi)*cos(chi)])


def stokes_vector(psi, chi, p=1, I=1):
    u"""A Stokes vector corresponding to a polarization ellipse with `psi`
    tilt, and `chi` circularity.

    INPUT:

    -  ``psi`` - The tilt of the polarization relative to the `x` axis.
    -  ``chi`` - The angle adjacent to the mayor axis of the polarization
        ellipse.
    -  ``p`` - The degree of polarization.
    -  ``I`` - The intensity of the field.

    OUTPUT:

    A Stokes vector.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> psi, chi, p, I = symbols("psi, chi, p, I", real=True)
    >>> pprint(stokes_vector(psi, chi, p, I), use_unicode=True)
    ⎡          I          ⎤
    ⎢                     ⎥
    ⎢I⋅p⋅cos(2⋅χ)⋅cos(2⋅ψ)⎥
    ⎢                     ⎥
    ⎢I⋅p⋅sin(2⋅ψ)⋅cos(2⋅χ)⎥
    ⎢                     ⎥
    ⎣    I⋅p⋅sin(2⋅χ)     ⎦


    Horizontal polarization
    >>> pprint(stokes_vector(0, 0), use_unicode=True)
    ⎡1⎤
    ⎢ ⎥
    ⎢1⎥
    ⎢ ⎥
    ⎢0⎥
    ⎢ ⎥
    ⎣0⎦

    Vertical polarization
    >>> pprint(stokes_vector(pi/2, 0), use_unicode=True)
    ⎡1 ⎤
    ⎢  ⎥
    ⎢-1⎥
    ⎢  ⎥
    ⎢0 ⎥
    ⎢  ⎥
    ⎣0 ⎦

    Diagonal polarization
    >>> pprint(stokes_vector(pi/4, 0), use_unicode=True)
    ⎡1⎤
    ⎢ ⎥
    ⎢0⎥
    ⎢ ⎥
    ⎢1⎥
    ⎢ ⎥
    ⎣0⎦

    Anti-diagonal polarization
    >>> pprint(stokes_vector(-pi/4, 0), use_unicode=True)
    ⎡1 ⎤
    ⎢  ⎥
    ⎢0 ⎥
    ⎢  ⎥
    ⎢-1⎥
    ⎢  ⎥
    ⎣0 ⎦

    Right-hand circular polarization
    >>> pprint(stokes_vector(0, pi/4), use_unicode=True)
    ⎡1⎤
    ⎢ ⎥
    ⎢0⎥
    ⎢ ⎥
    ⎢0⎥
    ⎢ ⎥
    ⎣1⎦

    Left-hand circular polarization
    >>> pprint(stokes_vector(0, -pi/4), use_unicode=True)
    ⎡1 ⎤
    ⎢  ⎥
    ⎢0 ⎥
    ⎢  ⎥
    ⎢0 ⎥
    ⎢  ⎥
    ⎣-1⎦

    Unpolarized light
    >>> pprint(stokes_vector(0, 0, 0), use_unicode=True)
    ⎡1⎤
    ⎢ ⎥
    ⎢0⎥
    ⎢ ⎥
    ⎢0⎥
    ⎢ ⎥
    ⎣0⎦

    """
    S0 = I
    S1 = I*p*cos(2*psi)*cos(2*chi)
    S2 = I*p*sin(2*psi)*cos(2*chi)
    S3 = I*p*sin(2*chi)
    return Matrix([S0, S1, S2, S3])


def jones_2_stokes(e):
    u"""Return the Stokes vector for a Jones vector `e`.

    INPUT:

    -  ``e`` - A Jones vector.

    OUTPUT:

    A Stokes vector.

    Examples
    ========

    >>> from sympy import pprint, Matrix
    >>> H = jones_vector(0, 0)
    >>> V = jones_vector(pi/2, 0)
    >>> D = jones_vector(pi/4, 0)
    >>> A = jones_vector(-pi/4, 0)
    >>> R = jones_vector(0, pi/4)
    >>> L = jones_vector(0, -pi/4)
    >>> pprint([jones_2_stokes(e) for e in [H, V, D, A, R, L]],
    ...         use_unicode=True)
    ⎡⎡1⎤  ⎡1 ⎤  ⎡1⎤  ⎡1 ⎤  ⎡1⎤  ⎡1 ⎤⎤
    ⎢⎢ ⎥  ⎢  ⎥  ⎢ ⎥  ⎢  ⎥  ⎢ ⎥  ⎢  ⎥⎥
    ⎢⎢1⎥  ⎢-1⎥  ⎢0⎥  ⎢0 ⎥  ⎢0⎥  ⎢0 ⎥⎥
    ⎢⎢ ⎥, ⎢  ⎥, ⎢ ⎥, ⎢  ⎥, ⎢ ⎥, ⎢  ⎥⎥
    ⎢⎢0⎥  ⎢0 ⎥  ⎢1⎥  ⎢-1⎥  ⎢0⎥  ⎢0 ⎥⎥
    ⎢⎢ ⎥  ⎢  ⎥  ⎢ ⎥  ⎢  ⎥  ⎢ ⎥  ⎢  ⎥⎥
    ⎣⎣0⎦  ⎣0 ⎦  ⎣0⎦  ⎣0 ⎦  ⎣1⎦  ⎣-1⎦⎦

    """
    ex, ey = e
    return Matrix([Abs(ex)**2 + Abs(ey)**2,
                   Abs(ex)**2 - Abs(ey)**2,
                   2*re(ex*ey.conjugate()),
                   -2*im(ex*ey.conjugate())])


def linear_polarizer(theta=0):
    u"""A linear polarizer Jones matrix with transmission axis at
    an angle `theta`.

    INPUT:

    -  ``theta`` - The angle of the transmission axis relative to the
        horizontal plane.

    OUTPUT:

    A Jones matrix representing the polarizer.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> theta = symbols("theta", real=True)
    >>> J = linear_polarizer(theta)
    >>> pprint(J)
    ⎡      2                     ⎤
    ⎢   cos (θ)     sin(θ)⋅cos(θ)⎥
    ⎢                            ⎥
    ⎢                     2      ⎥
    ⎣sin(θ)⋅cos(θ)     sin (θ)   ⎦


    """
    M = Matrix([[cos(theta)**2, sin(theta)*cos(theta)],
                [sin(theta)*cos(theta), sin(theta)**2]])
    return M


def phase_retarder(theta=0, delta=0):
    u"""A phase retarder Jones matrix with retardance `delta` at angle `theta`.

    INPUT:

    -  ``theta`` - The angle of the fast axis relative to the
        horizontal plane.
    -  ``delta`` - The phase difference between the fast and slow axes of the
        transmitted light.


    OUTPUT:

    A Jones matrix representing the retarder.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> theta, delta = symbols("theta, delta", real=True)
    >>> R = phase_retarder(theta, delta)
    >>> pprint(R, use_unicode=True)
    ⎡                          -ⅈ⋅δ               -ⅈ⋅δ               ⎤
    ⎢                          ─────              ─────              ⎥
    ⎢⎛ ⅈ⋅δ    2         2   ⎞    2    ⎛     ⅈ⋅δ⎞    2                ⎥
    ⎢⎝ℯ   ⋅sin (θ) + cos (θ)⎠⋅ℯ       ⎝1 - ℯ   ⎠⋅ℯ     ⋅sin(θ)⋅cos(θ)⎥
    ⎢                                                                ⎥
    ⎢            -ⅈ⋅δ                                           -ⅈ⋅δ ⎥
    ⎢            ─────                                          ─────⎥
    ⎢⎛     ⅈ⋅δ⎞    2                  ⎛ ⅈ⋅δ    2         2   ⎞    2  ⎥
    ⎣⎝1 - ℯ   ⎠⋅ℯ     ⋅sin(θ)⋅cos(θ)  ⎝ℯ   ⋅cos (θ) + sin (θ)⎠⋅ℯ     ⎦

    """
    R = Matrix([[cos(theta)**2 + exp(I*delta)*sin(theta)**2,
                (1-exp(I*delta))*cos(theta)*sin(theta)],
                [(1-exp(I*delta))*cos(theta)*sin(theta),
                sin(theta)**2 + exp(I*delta)*cos(theta)**2]])
    return R*exp(-I*delta/2)


def half_wave_retarder(theta):
    u"""A half-wave retarder Jones matrix at angle `theta`.

    INPUT:

    -  ``theta`` - The angle of the fast axis relative to the
        horizontal plane.


    OUTPUT:

    A Jones matrix representing the retarder.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> theta= symbols("theta", real=True)
    >>> HWP = half_wave_retarder(theta)
    >>> pprint(HWP, use_unicode=True)
    ⎡   ⎛     2         2   ⎞                        ⎤
    ⎢-ⅈ⋅⎝- sin (θ) + cos (θ)⎠    -2⋅ⅈ⋅sin(θ)⋅cos(θ)  ⎥
    ⎢                                                ⎥
    ⎢                             ⎛   2         2   ⎞⎥
    ⎣   -2⋅ⅈ⋅sin(θ)⋅cos(θ)     -ⅈ⋅⎝sin (θ) - cos (θ)⎠⎦

    """
    return phase_retarder(theta, pi)


def quarter_wave_retarder(theta):
    u"""A quarter-wave retarder Jones matrix at angle `theta`.

    INPUT:

    -  ``theta`` - The angle of the fast axis relative to the
        horizontal plane.


    OUTPUT:

    A Jones matrix representing the retarder.

    Examples
    ========

    A Jones matrix representing the retarder.

    >>> from sympy import pprint, symbols
    >>> theta= symbols("theta", real=True)
    >>> QWP = quarter_wave_retarder(theta)
    >>> pprint(QWP, use_unicode=True)
    ⎡                       -ⅈ⋅π            -ⅈ⋅π               ⎤
    ⎢                       ─────           ─────              ⎥
    ⎢⎛     2         2   ⎞    4               4                ⎥
    ⎢⎝ⅈ⋅sin (θ) + cos (θ)⎠⋅ℯ       (1 - ⅈ)⋅ℯ     ⋅sin(θ)⋅cos(θ)⎥
    ⎢                                                          ⎥
    ⎢         -ⅈ⋅π                                        -ⅈ⋅π ⎥
    ⎢         ─────                                       ─────⎥
    ⎢           4                  ⎛   2           2   ⎞    4  ⎥
    ⎣(1 - ⅈ)⋅ℯ     ⋅sin(θ)⋅cos(θ)  ⎝sin (θ) + ⅈ⋅cos (θ)⎠⋅ℯ     ⎦

    """
    return phase_retarder(theta, pi/2)


def attenuator(T):
    u"""An attenuator Jones matrix with transmissivity `T`.

    INPUT:

    -  ``T`` - The transmissivity of the attenuator.

    OUTPUT:

    A Jones matrix representing the attenuator.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> T = symbols("T", real=True)
    >>> NDF = attenuator(T)
    >>> pprint(NDF, use_unicode=True)
    ⎡√T  0 ⎤
    ⎢      ⎥
    ⎣0   √T⎦

    """
    return Matrix([[sqrt(T), 0], [0, sqrt(T)]])


def mirror(R):
    u"""A mirror Jones matrix with reflectivty `R`.

    INPUT:

    -  ``R`` - The transmissivity of the attenuator.

    OUTPUT:

    A Jones matrix representing the mirror.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> R = symbols("R", real=True)
    >>> pprint(mirror(R), use_unicode=True)
    ⎡√R   0 ⎤
    ⎢       ⎥
    ⎣0   -√R⎦

    """
    return Matrix([[sqrt(R), 0], [0, -sqrt(R)]])


def mueller_matrix(J):
    u"""A Mueller matrix corresponding to Jones matrix `J`.

    INPUT:

    -  ``J`` - A Jones matrix.

    OUTPUT:

    The corresponding Mueller matrix.

    Examples
    ========

    A Jones matrix representing the attenuator.

    >>> from sympy import pprint, symbols, pi, simplify
    >>> theta = symbols("theta", real=True)

    A linear_polarizer
    >>> pprint(mueller_matrix(linear_polarizer(theta)), use_unicode=True)
    ⎡            cos(2⋅θ)      sin(2⋅θ)     ⎤
    ⎢  1/2       ────────      ────────    0⎥
    ⎢               2             2         ⎥
    ⎢                                       ⎥
    ⎢cos(2⋅θ)  cos(4⋅θ)   1    sin(4⋅θ)     ⎥
    ⎢────────  ──────── + ─    ────────    0⎥
    ⎢   2         4       4       4         ⎥
    ⎢                                       ⎥
    ⎢sin(2⋅θ)    sin(4⋅θ)    1   cos(4⋅θ)   ⎥
    ⎢────────    ────────    ─ - ────────  0⎥
    ⎢   2           4        4      4       ⎥
    ⎢                                       ⎥
    ⎣   0           0             0        0⎦

    A half-wave plate
    >>> pprint(mueller_matrix(half_wave_retarder(theta)), use_unicode=True)
    ⎡1              0                           0               0 ⎤
    ⎢                                                             ⎥
    ⎢        4           2                                        ⎥
    ⎢0  8⋅sin (θ) - 8⋅sin (θ) + 1           sin(4⋅θ)            0 ⎥
    ⎢                                                             ⎥
    ⎢                                     4           2           ⎥
    ⎢0          sin(4⋅θ)           - 8⋅sin (θ) + 8⋅sin (θ) - 1  0 ⎥
    ⎢                                                             ⎥
    ⎣0              0                           0               -1⎦

    A quarter-wave plate
    >>> pprint(mueller_matrix(quarter_wave_retarder(theta)), use_unicode=True)
    ⎡1       0             0            0    ⎤
    ⎢                                        ⎥
    ⎢   cos(4⋅θ)   1    sin(4⋅θ)             ⎥
    ⎢0  ──────── + ─    ────────    -sin(2⋅θ)⎥
    ⎢      2       2       2                 ⎥
    ⎢                                        ⎥
    ⎢     sin(4⋅θ)    1   cos(4⋅θ)           ⎥
    ⎢0    ────────    ─ - ────────  cos(2⋅θ) ⎥
    ⎢        2        2      2               ⎥
    ⎢                                        ⎥
    ⎣0    sin(2⋅θ)     -cos(2⋅θ)        0    ⎦

    """
    A = Matrix([[1, 0, 0, 1],
                [1, 0, 0, -1],
                [0, 1, 1, 0],
                [0, -I, I, 0]])

    return simplify(A*TensorProduct(J, J.conjugate())*A.inv())


def polarizing_beam_splitter(Tp=1, Rs=1, Ts=0, Rp=0, phia=0, phib=0):
    u"""A polarizing beam splitter Jones matrix at angle `theta`.

    INPUT:

    -  ``Tp`` - The transmissivity of the P-polarized component.
    -  ``Rs`` - The reflectivity of the S-polarized component.
    -  ``Ts`` - The transmissivity of the S-polarized component.
    -  ``Rp`` - The reflectivity of the P-polarized component.
    -  ``phia`` - The phase difference between transmitted and reflected
                component for output mode a.
    -  ``phib`` - The phase difference between transmitted and reflected
                component for output mode b.

    OUTPUT:

    A Jones matrix representing the retarder.

    Examples
    ========

    >>> from sympy import pprint, symbols
    >>> Ts, Rs, Tp, Rp = symbols(r"Ts, Rs, Tp, Rp", positive=True)
    >>> phia, phib = symbols("phi_a, phi_b", real=True)
    >>> PBS = polarizing_beam_splitter(Tp, Rs, Ts, Rp, phia, phib)
    >>> pprint(PBS, use_unicode=True)
    ⎡   ____                        ____                ⎤
    ⎢ ╲╱ Tp           0         ⅈ⋅╲╱ Rp         0       ⎥
    ⎢                                                   ⎥
    ⎢                ____                     ____  ⅈ⋅φₐ⎥
    ⎢   0          ╲╱ Ts           0      ⅈ⋅╲╱ Rs ⋅ℯ    ⎥
    ⎢                                                   ⎥
    ⎢    ____                      ____                 ⎥
    ⎢ⅈ⋅╲╱ Rp          0          ╲╱ Tp          0       ⎥
    ⎢                                                   ⎥
    ⎢              ____  ⅈ⋅φ_b                  ____    ⎥
    ⎣   0      ⅈ⋅╲╱ Rs ⋅ℯ          0          ╲╱ Ts     ⎦

    """
    PBS = Matrix([[sqrt(Tp), 0, I*sqrt(Rp), 0],
                  [0, sqrt(Ts), 0, I*sqrt(Rs)*exp(I*phia)],
                  [I*sqrt(Rp), 0, sqrt(Tp), 0],
                  [0, I*sqrt(Rs)*exp(I*phib), 0, sqrt(Ts)]])
    return PBS
