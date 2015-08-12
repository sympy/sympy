# sympy/galgebra/tests/test_ga.py

"""
The reference D&L is "Geometric Algebra for Physicists" by Doran and Lasenby
"""

from sympy.core import expand, Rational, S, Symbol, symbols
from sympy.core.compatibility import range
from sympy.functions import sin, cos
from sympy.galgebra.ga import MV, Nga, Com
from sympy.galgebra.printing import GA_Printer
from sympy.matrices import Matrix
from sympy.simplify import collect, simplify
from sympy.utilities.pytest import XFAIL, slow


def F(x, n, nbar):
    """
    Conformal Mapping Function from 3D Euclidean space to 5D conformal space
    where the images of all maps are null vectors.
    """
    return Rational(1, 2)*((x*x)*n + 2*x - nbar)


def make_vector(a, m=3):
    global n, nbar
    if isinstance(a, str):
        sym_str = ''
        for i in range(m):
            sym_str += a + str(i + 1) + ' '
        sym_lst = list(symbols(sym_str))
        sym_lst.append(S.Zero)
        sym_lst.append(S.Zero)
        a = MV(sym_lst, 'vector')
    return F(a, n, nbar)


def test_rmul():
    """
    Test for commutative scalar multiplication.  Leftover from when sympy and
    numpy were not working together and __mul__ and __rmul__ would not give the
    same answer.
    """
    x, y, z = MV.setup('x y z')
    a, b, c = symbols('a b c')
    assert 5*x == x*5
    assert Rational(1, 2)*x == x*Rational(1, 2)
    assert a*x == x*a


def test_contraction():
    """
    Test for inner product and left and right contraction
    """

    e_1, e_2, e_3 = MV.setup('e_1 e_2 e_3', '1 0 0, 0 1 0, 0 0 1')

    assert ((e_1 ^ e_3) | e_1) == -e_3
    assert ((e_1 ^ e_3) > e_1) == -e_3
    assert (e_1 | (e_1 ^ e_3)) == e_3
    assert (e_1 < (e_1 ^ e_3)) == e_3
    assert ((e_1 ^ e_3) < e_1) == 0
    assert (e_1 > (e_1 ^ e_3)) == 0


def test_substitution():

    e_x, e_y, e_z = MV.setup('e_x e_y e_z', '1 0 0, 0 1 0, 0 0 1')
    x, y, z = symbols('x y z')

    X = x*e_x + y*e_y + z*e_z
    Y = X.subs([(x, 2), (y, 3), (z, 4)])
    assert Y == 2*e_x + 3*e_y + 4*e_z


def test_vector_extraction():
    """
    Show that conformal bivector encodes two points. See D&L Section 10.4.1
    """
    metric = ' 0 -1 #,' + \
             '-1  0 #,' + \
             ' #  # #,'

    P1, P2, a = MV.setup('P1 P2 a', metric)
    """
    P1 and P2 are null vectors and hence encode points in conformal space.
    Show that P1 and P2 can be extracted from the bivector B = P1^P2. a is a
    third vector in the conformal space with a.B not 0.
    """
    B = P1 ^ P2
    Bsq = B*B
    ap = a - (a ^ B)*B
    Ap = ap + ap*B
    Am = ap - ap*B
    P1dota = Symbol('(P1.a)')
    P2dota = Symbol('(P2.a)')
    Ap_test = (-2*P2dota)*P1
    Am_test = (-2*P1dota)*P2
    assert Ap == Ap_test
    assert Am == Am_test
    Ap2 = Ap*Ap
    Am2 = Am*Am
    assert Ap2 == S.Zero
    assert Am2 == S.Zero


def test_metrics():
    """
    Test specific metrics (diagpq, arbitrary_metric, arbitrary_metric_conformal)
    """
    from sympy.galgebra.ga import diagpq, arbitrary_metric
    metric = diagpq(3)
    p1, p2, p3 = MV.setup('p1 p2 p3', metric, debug=0)
    x1, y1, z1 = symbols('x1 y1 z1')
    x2, y2, z2 = symbols('x2 y2 z2')
    v1 = x1*p1 + y1*p2 + z1*p3
    v2 = x2*p1 + y2*p2 + z2*p3
    prod1 = v1*v2
    prod2 = (v1|v2) + (v1^v2)
    diff = prod1 - prod2
    assert diff == MV(S.Zero)
    metric = arbitrary_metric(3)
    p1, p2, p3 = MV.setup('p1 p2 p3', metric, debug=0)
    v1 = x1*p1 + y1*p2 + z1*p3
    v2 = x2*p1 + y2*p2 + z2*p3
    prod1 = v1*v2
    prod2 = (v1|v2) + (v1^v2)
    diff = prod1 - prod2
    assert diff == MV(S.Zero)


@XFAIL
def test_metrics_xfail():
    from sympy.galgebra.ga import arbitrary_metric_conformal
    metric = arbitrary_metric_conformal(3)
    p1, p2, p3 = MV.setup('p1 p2 p3', metric, debug=0)
    v1 = x1*p1 + y1*p2 + z1*p3
    v2 = x2*p1 + y2*p2 + z2*p3
    prod1 = v1*v2
    prod2 = (v1|v2) + (v1^v2)
    diff = prod1 - prod2
    assert diff == MV(S.Zero)


def test_geometry():
    """
    Test conformal geometric description of circles, lines, spheres, and planes.
    """
    metric = '1 0 0 0 0,' + \
             '0 1 0 0 0,' + \
             '0 0 1 0 0,' + \
             '0 0 0 0 2,' + \
             '0 0 0 2 0'

    e0, e1, e2, n, nbar = MV.setup('e0 e1 e2 n nbar', metric, debug=0)
    e = n + nbar
    #conformal representation of points

    A = F(e0, n, nbar)     # point a = (1,0,0)  A = F(a)
    B = F(e1, n, nbar)     # point b = (0,1,0)  B = F(b)
    C = F(-1*e0, n, nbar)  # point c = (-1,0,0) C = F(c)
    D = F(e2, n, nbar)     # point d = (0,0,1)  D = F(d)
    x0, x1, x2 = symbols('x0 x1 x2')
    X = F(MV([x0, x1, x2], 'vector'), n, nbar)

    Circle = A ^ B ^ C ^ X
    Line = A ^ B ^ n ^ X
    Sphere = A ^ B ^ C ^ D ^ X
    Plane = A ^ B ^ n ^ D ^ X

    #Circle through a, b, and c
    Circle_test = -x2*(e0 ^ e1 ^ e2 ^ n) + x2*(
        e0 ^ e1 ^ e2 ^ nbar) + Rational(1, 2)*(-1 + x0**2 + x1**2 + x2**2)*(e0 ^ e1 ^ n ^ nbar)
    diff = Circle - Circle_test
    assert diff == S.Zero

    #Line through a and b
    Line_test = -x2*(e0 ^ e1 ^ e2 ^ n) + \
        Rational(1, 2)*(-1 + x0 + x1)*(e0 ^ e1 ^ n ^ nbar) + \
        (Rational(1, 2)*x2)*(e0 ^ e2 ^ n ^ nbar) + \
        (-Rational(1, 2)*x2)*(e1 ^ e2 ^ n ^ nbar)
    diff = Line - Line_test
    assert diff == S.Zero

    #Sphere through a, b, c, and d
    Sphere_test = Rational(1, 2)*(1 - x0**2 - x1**2 - x2**2)*(e0 ^ e1 ^ e2 ^ n ^ nbar)
    diff = Sphere - Sphere_test
    assert diff == S.Zero

    #Plane through a, b, and d
    Plane_test = Rational(1, 2)*(1 - x0 - x1 - x2)*(e0 ^ e1 ^ e2 ^ n ^ nbar)
    diff = Plane - Plane_test
    assert diff == S.Zero


@slow
def test_extract_plane_and_line():
    """
    Show that conformal trivector encodes planes and lines. See D&L section
    10.4.2
    """
    metric = '# # # 0 0,' + \
             '# # # 0 0,' + \
             '# # # 0 0,' + \
             '0 0 0 0 2,' + \
             '0 0 0 2 0'

    p1, p2, p3, n, nbar = MV.setup('p1 p2 p3 n nbar', metric, debug=0)

    P1 = F(p1, n, nbar)
    P2 = F(p2, n, nbar)
    P3 = F(p3, n, nbar)

    #Line through p1 and p2
    L = P1 ^ P2 ^ n
    delta = (L | n) | nbar
    delta_test = 2*p1 - 2*p2
    diff = delta - delta_test
    assert diff == S.Zero

    #Plane through p1, p2, and p3
    C = P1 ^ P2 ^ P3
    delta = ((C ^ n) | n) | nbar
    delta_test = 2*(p1 ^ p2) - 2*(p1 ^ p3) + 2*(p2 ^ p3)
    diff = delta - delta_test
    assert diff == S.Zero


@XFAIL
def test_reciprocal_frame():
    """
    Test of formula for general reciprocal frame of three vectors.
    Let three independent vectors be e1, e2, and e3. The reciprocal
    vectors E1, E2, and E3 obey the relations:

    e_i.E_j = delta_ij*(e1^e2^e3)**2
    """
    metric = '1 # #,' + \
             '# 1 #,' + \
             '# # 1,'

    e1, e2, e3 = MV.setup('e1 e2 e3', metric)
    E = e1 ^ e2 ^ e3
    Esq = (E*E)()
    Esq_inv = 1/Esq
    E1 = (e2 ^ e3)*E
    E2 = (-1)*(e1 ^ e3)*E
    E3 = (e1 ^ e2)*E
    w = (E1 | e2)
    w.collect(MV.g)
    w = w().expand()
    w = (E1 | e3)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E2 | e1)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E2 | e3)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E3 | e1)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E3 | e2)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E1 | e1)
    w = w().expand()
    Esq = Esq.expand()
    assert w/Esq == 1
    w = (E2 | e2)
    w = w().expand()
    assert w/Esq == 1
    w = (E3 | e3)
    w = w().expand()
    assert w/Esq == 1


@XFAIL
def test_derivative():
    coords = x, y, z = symbols('x y z')
    e_x, e_y, e_z, _ = MV.setup('e', '1 0 0, 0 1 0, 0 0 1', coords=coords)
    X = x*e_x + y*e_y + z*e_z
    a = MV('a', 'vector')

    assert ((X | a).grad()) == a
    assert ((X*X).grad()) == 2*X
    assert (X*X*X).grad() == 5*X*X
    assert X.grad_int() == 3


@XFAIL
def test_str():
    e_1, e_2, e_3 = MV.setup('e_1 e_2 e_3', '1 0 0, 0 1 0, 0 0 1')

    X = MV('x')
    assert str(X) == 'x + x__1*e_1 + x__2*e_2 + x__3*e_3 + x__12*e_1^e_2 + x__13*e_1^e_3 + x__23*e_2^e_3 + x__123**e_1^e_2^e_3'
    Y = MV('y', 'spinor')
    assert str(Y) == 'y + y__12*e_1^e_2 + y__13*e_1^e_3 + y__23*e_2^e_3'
    Z = X + Y
    assert str(Z) == 'x + y + x__1*e_1 + x__2*e_2 + x__3*e_3 + (x__12 + y__12)*e_1^e_2 + (x__13 + y__13)*e_1^e_3 + (x__23 + y__23)*e_2^e_3 + x__123*e_1^e_2^e_3'
    assert str(e_1 | e_1) == '1'


@XFAIL
def test_metric():
    MV.setup('e_1 e_2 e_3', '[1,1,1]')
    assert MV.metric == Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


@XFAIL
def test_constructor():
    """
    Test various multivector constructors
    """
    e_1, e_2, e_3 = MV.setup('e_1 e_2 e_3', '[1,1,1]')
    assert str(MV('a', 'scalar')) == 'a'
    assert str(MV('a', 'vector')) == 'a__1*e_1 + a__2*e_2 + a__3*e_3'
    assert str(MV('a', 'pseudo')) == 'a__123*e_1^e_2^e_3'
    assert str(MV('a', 'spinor')) == 'a + a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3'
    assert str(MV('a')) == 'a + a__1*e_1 + a__2*e_2 + a__3*e_3 + a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3 + a__123*e_1^e_2^e_3'
    assert str(MV([2, 'a'], 'grade')) == 'a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3'
    assert str(MV('a', 'grade2')) == 'a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3'


def test_basic_multivector_operations():
    with GA_Printer():
        (ex, ey, ez) = MV.setup('e*x|y|z')

        A = MV('A', 'mv')

        assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'
        assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'
        assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'

        X = MV('X', 'vector')
        Y = MV('Y', 'vector')

        assert str(X) == 'X__x*e_x + X__y*e_y + X__z*e_z'
        assert str(Y) == 'Y__x*e_x + Y__y*e_y + Y__z*e_z'

        assert str((X*Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z + (X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
        assert str((X ^ Y)) == '(X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
        assert str((X | Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z'

        (ex, ey) = MV.setup('e*x|y')

        X = MV('X', 'vector')
        A = MV('A', 'spinor')

        assert str(X) == 'X__x*e_x + X__y*e_y'
        assert str(A) == 'A + A__xy*e_x^e_y'

        assert str((X | A)) == '(-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y))*e_x + (A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y))*e_y'
        assert str((X < A)) == '(-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y))*e_x + (A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y))*e_y'
        assert str((A > X)) == '(A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y))*e_x + (-A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y))*e_y'

        (ex, ey) = MV.setup('e*x|y', metric='[1,1]')

        X = MV('X', 'vector')
        A = MV('A', 'spinor')

        assert str(X) == 'X__x*e_x + X__y*e_y'
        assert str(A) == 'A + A__xy*e_x^e_y'

        assert str((X*A)) == '(A*X__x - A__xy*X__y)*e_x + (A*X__y + A__xy*X__x)*e_y'
        assert str((X | A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
        assert str((X < A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
        assert str((X > A)) == 'A*X__x*e_x + A*X__y*e_y'

        assert str((A*X)) == '(A*X__x + A__xy*X__y)*e_x + (A*X__y - A__xy*X__x)*e_y'
        assert str((A | X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'
        assert str((A < X)) == 'A*X__x*e_x + A*X__y*e_y'
        assert str((A > X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'

    return


@slow
def test_check_generalized_BAC_CAB_formulas():
    with GA_Printer():
        (a, b, c, d, e) = MV.setup('a b c d e')

        assert str(a | (b*c)) == '-(a.c)*b + (a.b)*c'
        assert str(a | (b ^ c)) == '-(a.c)*b + (a.b)*c'
        assert str(a | (b ^ c ^ d)) == '(a.d)*b^c - (a.c)*b^d + (a.b)*c^d'
        assert str((a | (b ^ c)) + (c | (a ^ b)) + (b | (c ^ a))) == '0'
        assert str(a*(b ^ c) - b*(a ^ c) + c*(a ^ b)) == '3*a^b^c'
        assert str(a*(b ^ c ^ d) - b*(a ^ c ^ d) + c*(a ^ b ^ d) - d*(a ^ b ^ c)) == '4*a^b^c^d'
        assert str((a ^ b) | (c ^ d)) == '-(a.c)*(b.d) + (a.d)*(b.c)'
        assert str(((a ^ b) | c) | d) == '-(a.c)*(b.d) + (a.d)*(b.c)'
        assert str(Com(a ^ b, c ^ d)) == '-(b.d)*a^c + (b.c)*a^d + (a.d)*b^c - (a.c)*b^d'
        assert str((a | (b ^ c)) | (d ^ e)) == '(-(a.b)*(c.e) + (a.c)*(b.e))*d + ((a.b)*(c.d) - (a.c)*(b.d))*e'

    return

def test_derivatives_in_rectangular_coordinates():
    with GA_Printer():
        X = (x, y, z) = symbols('x y z')
        (ex, ey, ez, grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=X)

        f = MV('f', 'scalar', fct=True)
        A = MV('A', 'vector', fct=True)
        B = MV('B', 'grade2', fct=True)
        C = MV('C', 'mv', fct=True)

        assert str(f) == 'f'
        assert str(A) == 'A__x*e_x + A__y*e_y + A__z*e_z'
        assert str(B) == 'B__xy*e_x^e_y + B__xz*e_x^e_z + B__yz*e_y^e_z'
        assert str(C) == 'C + C__x*e_x + C__y*e_y + C__z*e_z + C__xy*e_x^e_y + C__xz*e_x^e_z + C__yz*e_y^e_z + C__xyz*e_x^e_y^e_z'

        assert str(grad*f) == 'D{x}f*e_x + D{y}f*e_y + D{z}f*e_z'
        assert str(grad | A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
        assert str(grad*A) == 'D{x}A__x + D{y}A__y + D{z}A__z + (-D{y}A__x + D{x}A__y)*e_x^e_y + (-D{z}A__x + D{x}A__z)*e_x^e_z + (-D{z}A__y + D{y}A__z)*e_y^e_z'

        assert str(-MV.I*(grad ^ A)) == '(-D{z}A__y + D{y}A__z)*e_x + (D{z}A__x - D{x}A__z)*e_y + (-D{y}A__x + D{x}A__y)*e_z'
        assert str(grad*B) == '(-(D{y}B__xy + D{z}B__xz))*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z + (D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
        assert str(grad ^ B) == '(D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
        assert str(grad | B) == '(-(D{y}B__xy + D{z}B__xz))*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'

        assert str(grad < A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
        assert str(grad > A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
        assert str(grad < B) == '(-(D{y}B__xy + D{z}B__xz))*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'
        assert str(grad > B) == '0'
        assert str(grad < C) == 'D{x}C__x + D{y}C__y + D{z}C__z + (-(D{y}C__xy + D{z}C__xz))*e_x + (D{x}C__xy - D{z}C__yz)*e_y + (D{x}C__xz + D{y}C__yz)*e_z + D{z}C__xyz*e_x^e_y - D{y}C__xyz*e_x^e_z + D{x}C__xyz*e_y^e_z'
        assert str(grad > C) == 'D{x}C__x + D{y}C__y + D{z}C__z + D{x}C*e_x + D{y}C*e_y + D{z}C*e_z'

    return

def test_derivatives_in_spherical_coordinates():
    with GA_Printer():
        X = (r, th, phi) = symbols('r theta phi')
        curv = [[r*cos(phi)*sin(th), r*sin(phi)*sin(th), r*cos(th)], [1, r, r*sin(th)]]
        (er, eth, ephi, grad) = MV.setup('e_r e_theta e_phi', metric='[1,1,1]', coords=X, curv=curv)

        f = MV('f', 'scalar', fct=True)
        A = MV('A', 'vector', fct=True)
        B = MV('B', 'grade2', fct=True)

        assert str(f) == 'f'
        assert str(A) == 'A__r*e_r + A__theta*e_theta + A__phi*e_phi'
        assert str(B) == 'B__rtheta*e_r^e_theta + B__rphi*e_r^e_phi + B__thetaphi*e_theta^e_phi'

        assert str(grad*f) == 'D{r}f*e_r + D{theta}f/r*e_theta + D{phi}f/(r*sin(theta))*e_phi'
        assert str(grad | A) == 'D{r}A__r + 2*A__r/r + A__theta*cos(theta)/(r*sin(theta)) + D{theta}A__theta/r + D{phi}A__phi/(r*sin(theta))'
        assert str(-MV.I*(grad ^ A)) == '((A__phi*cos(theta)/sin(theta) + D{theta}A__phi - D{phi}A__theta/sin(theta))/r)*e_r + (-D{r}A__phi - A__phi/r + D{phi}A__r/(r*sin(theta)))*e_theta + (D{r}A__theta + A__theta/r - D{theta}A__r/r)*e_phi'
        assert str(grad ^ B) == '(D{r}B__thetaphi - B__rphi*cos(theta)/(r*sin(theta)) + 2*B__thetaphi/r - D{theta}B__rphi/r + D{phi}B__rtheta/(r*sin(theta)))*e_r^e_theta^e_phi'

    return

def test_rounding_numerical_components():
    with GA_Printer():
        (ex, ey, ez) = MV.setup('e_x e_y e_z', metric='[1,1,1]')

        X = 1.2*ex + 2.34*ey + 0.555*ez
        Y = 0.333*ex + 4*ey + 5.3*ez

        assert str(X) == '1.20000000000000*e_x + 2.34000000000000*e_y + 0.555000000000000*e_z'
        assert str(Nga(X, 2)) == '1.2*e_x + 2.3*e_y + 0.55*e_z'
        assert str(X*Y) == '12.7011000000000 + 4.02078000000000*e_x^e_y + 6.17518500000000*e_x^e_z + 10.1820000000000*e_y^e_z'
        assert str(Nga(X*Y, 2)) == '13. + 4.0*e_x^e_y + 6.2*e_x^e_z + 10.*e_y^e_z'

    return

def test_noneuclidian_distance_calculation():
    from sympy import solve, sqrt
    from sympy.solvers.solveset import solveset
    with GA_Printer():
        metric = '0 # #,# 0 #,# # 1'
        (X, Y, e) = MV.setup('X Y e', metric)

        assert str((X ^ Y)*(X ^ Y)) == '(X.Y)**2'

        L = X ^ Y ^ e
        B = L*e
        assert str(B) == 'X^Y - (Y.e)*X^e + (X.e)*Y^e'
        Bsq = B*B
        assert str(Bsq) == '(X.Y)*((X.Y) - 2*(X.e)*(Y.e))'
        Bsq = Bsq.scalar()
        assert str(B) == 'X^Y - (Y.e)*X^e + (X.e)*Y^e'

        BeBr = B*e*B.rev()
        assert str(BeBr) == '((X.Y)*(-(X.Y) + 2*(X.e)*(Y.e)))*e'
        assert str(B*B) == '(X.Y)*((X.Y) - 2*(X.e)*(Y.e))'
        assert str(L*L) == '(X.Y)*((X.Y) - 2*(X.e)*(Y.e))'
        (s, c, Binv, M, BigS, BigC, alpha, XdotY, Xdote, Ydote) = symbols('s c (1/B) M S C alpha (X.Y) (X.e) (Y.e)')

        Bhat = Binv*B
        R = c + s*Bhat
        assert str(R) == 'c + (1/B)*s*X^Y - (1/B)*(Y.e)*s*X^e + (1/B)*(X.e)*s*Y^e'

        Z = R*X*R.rev()
        Z.obj = expand(Z.obj)
        Z.obj = Z.obj.collect([Binv, s, c, XdotY])
        assert str(Z) == '((1/B)**2*(X.Y)**2*s**2 - 2*(1/B)**2*(X.Y)*(X.e)*(Y.e)*s**2 + 2*(1/B)*(X.Y)*c*s - 2*(1/B)*(X.e)*(Y.e)*c*s + c**2)*X + 2*(1/B)*(X.e)**2*c*s*Y + (2*(1/B)*(X.Y)*(X.e)*s*(-(1/B)*(X.Y)*s + 2*(1/B)*(X.e)*(Y.e)*s - c))*e'
        W = Z | Y
        # From this point forward all calculations are with sympy scalars
        W = W.scalar()
        assert str(W) == '(1/B)**2*(X.Y)**3*s**2 - 4*(1/B)**2*(X.Y)**2*(X.e)*(Y.e)*s**2 + 4*(1/B)**2*(X.Y)*(X.e)**2*(Y.e)**2*s**2 + 2*(1/B)*(X.Y)**2*c*s - 4*(1/B)*(X.Y)*(X.e)*(Y.e)*c*s + (X.Y)*c**2'
        W = expand(W)
        W = simplify(W)
        W = W.collect([s*Binv])

        M = 1/Bsq
        W = W.subs(Binv**2, M)
        W = simplify(W)
        Bmag = sqrt(XdotY**2 - 2*XdotY*Xdote*Ydote)
        W = W.collect([Binv*c*s, XdotY])

        #Double angle substitutions

        W = W.subs(2*XdotY**2 - 4*XdotY*Xdote*Ydote, 2/(Binv**2))
        W = W.subs(2*c*s, BigS)
        W = W.subs(c**2, (BigC + 1)/2)
        W = W.subs(s**2, (BigC - 1)/2)
        W = simplify(W)
        W = expand(W)
        W = W.subs(1/Binv, Bmag)

        assert str(W) == '(X.Y)*C - (X.e)*(Y.e)*C + (X.e)*(Y.e) + S*sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'

        Wd = collect(W, [BigC, BigS], exact=True, evaluate=False)

        Wd_1 = Wd[S.One]
        Wd_C = Wd[BigC]
        Wd_S = Wd[BigS]

        assert str(Wd_1) == '(X.e)*(Y.e)'
        assert str(Wd_C) == '(X.Y) - (X.e)*(Y.e)'
        assert str(Wd_S) == 'sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'

        assert str(Bmag) == 'sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'
        Wd_1 = Wd_1.subs(Bmag, 1/Binv)
        Wd_C = Wd_C.subs(Bmag, 1/Binv)
        Wd_S = Wd_S.subs(Bmag, 1/Binv)

        lhs = Wd_1 + Wd_C*BigC
        rhs = -Wd_S*BigS
        lhs = lhs**2
        rhs = rhs**2
        W = expand(lhs - rhs)
        W = expand(W.subs(1/Binv**2, Bmag**2))
        W = expand(W.subs(BigS**2, BigC**2 - 1))
        W = W.collect([BigC, BigC**2], evaluate=False)

        a = simplify(W[BigC**2])
        b = simplify(W[BigC])
        c = simplify(W[S.One])

        assert str(a) == '(X.e)**2*(Y.e)**2'
        assert str(b) == '2*(X.e)*(Y.e)*((X.Y) - (X.e)*(Y.e))'
        assert str(c) == '(X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e) + (X.e)**2*(Y.e)**2'

        x = Symbol('x')
        C1 = solve(a*x**2 + b*x + c, x)[0]
        C2 = solveset(a*x**2 + b*x + c, x).args[0]
        assert str(expand(simplify(expand(C1)))) == '-(X.Y)/((X.e)*(Y.e)) + 1'
        assert str(expand(simplify(expand(C2)))) == '-(X.Y)/((X.e)*(Y.e)) + 1'

    return

def test_conformal_representations_of_circles_lines_spheres_and_planes():
    global n, nbar
    with GA_Printer():

        metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

        (e1, e2, e3, n, nbar) = MV.setup('e_1 e_2 e_3 n nbar', metric)

        e = n + nbar
        #conformal representation of points

        A = make_vector(e1)
        B = make_vector(e2)
        C = make_vector(-e1)
        D = make_vector(e3)
        X = make_vector('x', 3)

        assert str(A) == 'e_1 + 1/2*n - 1/2*nbar'
        assert str(B) == 'e_2 + 1/2*n - 1/2*nbar'
        assert str(C) == '-e_1 + 1/2*n - 1/2*nbar'
        assert str(D) == 'e_3 + 1/2*n - 1/2*nbar'
        assert str(X) == 'x1*e_1 + x2*e_2 + x3*e_3 + ((x1**2 + x2**2 + x3**2)/2)*n - 1/2*nbar'

        assert str((A ^ B ^ C ^ X)) == '-x3*e_1^e_2^e_3^n + x3*e_1^e_2^e_3^nbar + ((x1**2 + x2**2 + x3**2 - 1)/2)*e_1^e_2^n^nbar'
        assert str((A ^ B ^ n ^ X)) == '-x3*e_1^e_2^e_3^n + ((x1 + x2 - 1)/2)*e_1^e_2^n^nbar + x3/2*e_1^e_3^n^nbar - x3/2*e_2^e_3^n^nbar'
        assert str((((A ^ B) ^ C) ^ D) ^ X) == '((-x1**2 - x2**2 - x3**2 + 1)/2)*e_1^e_2^e_3^n^nbar'
        assert str((A ^ B ^ n ^ D ^ X)) == '((-x1 - x2 - x3 + 1)/2)*e_1^e_2^e_3^n^nbar'

        L = (A ^ B ^ e) ^ X

        assert str(L) == '-x3*e_1^e_2^e_3^n - x3*e_1^e_2^e_3^nbar + (-x1**2/2 + x1 - x2**2/2 + x2 - x3**2/2 - 1/2)*e_1^e_2^n^nbar + x3*e_1^e_3^n^nbar - x3*e_2^e_3^n^nbar'

    return


@slow
def test_properties_of_geometric_objects():
    with GA_Printer():
        metric = '# # # 0 0,' + \
                 '# # # 0 0,' + \
                 '# # # 0 0,' + \
                 '0 0 0 0 2,' + \
                 '0 0 0 2 0'

        (p1, p2, p3, n, nbar) = MV.setup('p1 p2 p3 n nbar', metric)

        P1 = F(p1, n, nbar)
        P2 = F(p2, n, nbar)
        P3 = F(p3, n, nbar)

        L = P1 ^ P2 ^ n
        delta = (L | n) | nbar
        assert str(delta) == '2*p1 - 2*p2'

        C = P1 ^ P2 ^ P3
        delta = ((C ^ n) | n) | nbar
        assert str(delta) == '2*p1^p2 - 2*p1^p3 + 2*p2^p3'
        assert str((p2 - p1) ^ (p3 - p1)) == 'p1^p2 - p1^p3 + p2^p3'

    return

def test_extracting_vectors_from_conformal_2_blade():
    with GA_Printer():
        metric = ' 0 -1 #,' + \
                 '-1  0 #,' + \
                 ' #  # #,'

        (P1, P2, a) = MV.setup('P1 P2 a', metric)

        B = P1 ^ P2
        Bsq = B*B
        assert str(Bsq) == '1'
        ap = a - (a ^ B)*B
        assert str(ap) == '-(P2.a)*P1 - (P1.a)*P2'

        Ap = ap + ap*B
        Am = ap - ap*B

        assert str(Ap) == '-2*(P2.a)*P1'
        assert str(Am) == '-2*(P1.a)*P2'

        assert str(Ap*Ap) == '0'
        assert str(Am*Am) == '0'

        aB = a | B
        assert str(aB) == '-(P2.a)*P1 + (P1.a)*P2'

    return

def test_reciprocal_frame_test():
    with GA_Printer():
        metric = '1 # #,' + \
                 '# 1 #,' + \
                 '# # 1,'

        (e1, e2, e3) = MV.setup('e1 e2 e3', metric)

        E = e1 ^ e2 ^ e3
        Esq = (E*E).scalar()
        assert str(E) == 'e1^e2^e3'
        assert str(Esq) == '(e1.e2)**2 - 2*(e1.e2)*(e1.e3)*(e2.e3) + (e1.e3)**2 + (e2.e3)**2 - 1'
        Esq_inv = 1/Esq

        E1 = (e2 ^ e3)*E
        E2 = (-1)*(e1 ^ e3)*E
        E3 = (e1 ^ e2)*E

        assert str(E1) == '((e2.e3)**2 - 1)*e1 + ((e1.e2) - (e1.e3)*(e2.e3))*e2 + (-(e1.e2)*(e2.e3) + (e1.e3))*e3'
        assert str(E2) == '((e1.e2) - (e1.e3)*(e2.e3))*e1 + ((e1.e3)**2 - 1)*e2 + (-(e1.e2)*(e1.e3) + (e2.e3))*e3'
        assert str(E3) == '(-(e1.e2)*(e2.e3) + (e1.e3))*e1 + (-(e1.e2)*(e1.e3) + (e2.e3))*e2 + ((e1.e2)**2 - 1)*e3'

        w = (E1 | e2)
        w = w.expand()
        assert str(w) == '0'

        w = (E1 | e3)
        w = w.expand()
        assert str(w) == '0'

        w = (E2 | e1)
        w = w.expand()
        assert str(w) == '0'

        w = (E2 | e3)
        w = w.expand()
        assert str(w) == '0'

        w = (E3 | e1)
        w = w.expand()
        assert str(w) == '0'

        w = (E3 | e2)
        w = w.expand()
        assert str(w) == '0'

        w = (E1 | e1)
        w = (w.expand()).scalar()
        Esq = expand(Esq)
        assert str(simplify(w/Esq)) == '1'

        w = (E2 | e2)
        w = (w.expand()).scalar()
        assert str(simplify(w/Esq)) == '1'

        w = (E3 | e3)
        w = (w.expand()).scalar()
        assert str(simplify(w/Esq)) == '1'

    return
