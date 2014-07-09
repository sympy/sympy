# sympy/galgebra/tests/test_ga.py

"""
The reference D&L is "Geometric Algebra for Physicists" by Doran and Lasenby
"""

from sympy.core import expand, Rational, S, Symbol, symbols
from sympy.functions import sin, cos
from sympy.galgebra.ga import Ga
from sympy.galgebra.mv import Com, Nga
from sympy.galgebra.lt import Mlt
from sympy.matrices import Matrix
from sympy.simplify import collect, simplify
from sympy.utilities.pytest import XFAIL
from sympy.galgebra.metric import linear_expand
from sympy import diff, Function, expand


def F(x, n, nbar):
    """
    Conformal Mapping Function from 3D Euclidean space to 5D conformal space
    where the images of all maps are null vectors.
    """
    return Rational(1, 2)*((x*x)*n + 2*x - nbar)

def make_vector(a, ga):
    global n, nbar
    if isinstance(a, str):
        a = ga.mv(a, 'vector')
        a.set_coef(1,3,S(0))
        a.set_coef(1,4,S(0))
    return F(a, n, nbar)

def test_differential_operators():

    xyz_coords = (x, y, z) = symbols('x y z', real=True)
    (o3d, ex, ey, ez) = Ga.build('e', g=[1, 1, 1], coords=xyz_coords)
    f = o3d.mv('f', 'scalar', f=True)
    lap = o3d.grad*o3d.grad

    assert str(lap) == 'D{x}^2 + D{y}^2 + D{z}^2'
    assert str(lap * f) == 'D{x}^2f + D{y}^2f + D{z}^2f'

    sph_coords = (r, th, phi) = symbols('r theta phi', real=True)
    (sp3d, er, eth, ephi) = Ga.build('e', g=[1, r**2, r**2 * sin(th)**2], coords=sph_coords, norm=True)
    f = sp3d.mv('f', 'scalar', f=True)
    lap = sp3d.grad*sp3d.grad
    assert str(lap) == '2/r*D{r} + cos(theta)/(r**2*sin(theta))*D{theta} + D{r}^2 + r**(-2)*D{theta}^2 + 1/(r**2*sin(theta)**2)*D{phi}^2'
    assert str(lap * f) == 'D{r}^2f + 2*D{r}f/r + D{theta}^2f/r**2 + cos(theta)*D{theta}f/(r**2*sin(theta)) + D{phi}^2f/(r**2*sin(theta)**2)'

    A = o3d.mv('A','vector')
    xs = o3d.mv(x)

    assert o3d.grad*A == 0
    assert str(A*o3d.grad) == 'A__x*D{x} + A__y*D{y} + A__z*D{z} + e_x^e_y*(-A__y*D{x} + A__x*D{y}) + e_x^e_z*(-A__z*D{x} + A__x*D{z}) + e_y^e_z*(-A__z*D{y} + A__y*D{z})'
    assert o3d.grad*xs == ex
    assert str(xs*o3d.grad) == 'e_x*x*D{x} + e_y*x*D{y} + e_z*x*D{z}'
    assert str(o3d.grad*(o3d.grad+xs)) == 'D{x}^2 + D{y}^2 + D{z}^2 + e_x*D{}'
    assert str((o3d.grad+xs)*o3d.grad) == 'D{x}^2 + D{y}^2 + D{z}^2 + e_x*x*D{x} + e_y*x*D{y} + e_z*x*D{z}'

    return

def test_rmul():
    """
    Test for commutative scalar multiplication.  Leftover from when sympy and
    numpy were not working together and __mul__ and __rmul__ would not give the
    same answer.
    """
    g3d, x, y, z = Ga.build('x y z')
    a, b, c = symbols('a b c')
    assert 5*x == x*5
    assert Rational(1, 2)*x == x*Rational(1, 2)
    assert a*x == x*a


def test_contraction():
    """
    Test for inner product and left and right contraction
    """

    o3d, e_1, e_2, e_3 = Ga.build('e_1 e_2 e_3', g='1 0 0, 0 1 0, 0 0 1')

    assert ((e_1 ^ e_3) | e_1) == -e_3
    assert ((e_1 ^ e_3) > e_1) == -e_3
    assert (e_1 | (e_1 ^ e_3)) == e_3
    assert (e_1 < (e_1 ^ e_3)) == e_3
    assert ((e_1 ^ e_3) < e_1) == 0
    assert (e_1 > (e_1 ^ e_3)) == 0


def test_substitution():

    o3d, e_x, e_y, e_z = Ga.build('e_x e_y e_z', g='1 0 0, 0 1 0, 0 0 1')
    x, y, z = symbols('x y z')

    X = x*e_x + y*e_y + z*e_z
    Y = X.subs([(x, 2), (y, 3), (z, 4)])
    assert Y == 2*e_x + 3*e_y + 4*e_z


def test_vector_extraction():
    """
    Show that conformal bivector encodes two points. See D&L Section 10.4.1
    """
    metric = '0 -1 #,' + \
             '-1 0 #,' + \
             '# # #'

    cext, P1, P2, a = Ga.build('P1 P2 a', g=metric)
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

    P1dota = cext.g[0,2]
    P2dota = cext.g[1,2]
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

    o3d, p1, p2, p3 = Ga.build('p1 p2 p3', g=[1,1,1], debug=0)
    x1, y1, z1 = symbols('x1 y1 z1')
    x2, y2, z2 = symbols('x2 y2 z2')
    v1 = x1*p1 + y1*p2 + z1*p3
    v2 = x2*p1 + y2*p2 + z2*p3
    prod1 = v1*v2
    prod2 = (v1|v2) + (v1^v2)
    diff = prod1 - prod2
    assert diff == o3d.mv(S.Zero)

    g3d, p1, p2, p3 = Ga.build('p1 p2 p3', g='# # #, # # #, # # #', debug=0)
    v1 = x1*p1 + y1*p2 + z1*p3
    v2 = x2*p1 + y2*p2 + z2*p3
    prod1 = v1*v2
    prod2 = (v1|v2) + (v1^v2)
    diff = prod1 - prod2
    assert diff == g3d.mv(S.Zero)


def test_geometry():
    """
    Test conformal geometric description of circles, lines, spheres, and planes.
    """
    metric = '1 0 0 0 0,' + \
             '0 1 0 0 0,' + \
             '0 0 1 0 0,' + \
             '0 0 0 0 2,' + \
             '0 0 0 2 0'

    cf3d, e0, e1, e2, n, nbar = Ga.build('e0 e1 e2 n nbar', g=metric)
    e = n + nbar
    #conformal representation of points

    A = F(e0, n, nbar)     # point a = (1,0,0)  A = F(a)
    B = F(e1, n, nbar)     # point b = (0,1,0)  B = F(b)
    C = F(-1*e0, n, nbar)  # point c = (-1,0,0) C = F(c)
    D = F(e2, n, nbar)     # point d = (0,0,1)  D = F(d)
    x0, x1, x2 = symbols('x0 x1 x2')
    x = x0 * e0 + x1 * e1 + x2 * e2
    X = F(x, n, nbar)

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

    cf3d, p1, p2, p3, n, nbar = Ga.build('p1 p2 p3 n nbar', g=metric)

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


def test_reciprocal_frame():
    """
    Test of formula for general reciprocal frame of three vectors.
    Let three independent vectors be e1, e2, and e3. The reciprocal
    vectors E1, E2, and E3 obey the relations:

    e_i.E_j = delta_ij*(e1^e2^e3)**2
    """
    g = '1 # #,'+ \
        '# 1 #,'+ \
        '# # 1'

    g3dn = Ga('e1 e2 e3',g=g)

    (e1,e2,e3) = g3dn.mv()

    E = e1^e2^e3
    Esq = (E*E).scalar()
    Esq_inv = 1 / Esq

    E1 = (e2^e3)*E
    E2 = (-1)*(e1^e3)*E
    E3 = (e1^e2)*E

    w = (E1|e2)
    w = w.expand()
    assert w.scalar() == 0

    w = (E1|e3)
    w = w.expand()
    assert w.scalar() == 0

    w = (E2|e1)
    w = w.expand()
    assert w.scalar() == 0

    w = (E2|e3)
    w = w.expand()
    assert w.scalar() == 0

    w = (E3|e1)
    w = w.expand()
    assert w.scalar() == 0

    w = (E3|e2)
    w = w.expand()
    assert w.scalar() == 0

    w = (E1|e1)
    w = (w.expand()).scalar()
    Esq = expand(Esq)
    assert simplify(w/Esq) == 1

    w = (E2|e2)
    w = (w.expand()).scalar()
    assert simplify(w/Esq) == 1

    w = (E3|e3)
    w = (w.expand()).scalar()
    assert simplify(w/Esq) == 1


def test_derivative():
    coords = x, y, z = symbols('x y z')
    o3d, e_x, e_y, e_z = Ga.build('e', g=[1, 1, 1], coords=coords)
    grad = o3d.grad
    X = x*e_x + y*e_y + z*e_z
    a = o3d.mv('a', 'vector')

    assert (grad * (X | a)) == a
    assert (grad * (X*X)) == 2*X
    assert grad * (X*X*X) == 5*X*X
    assert (grad *X).scalar() == 3


def test_metric():
    o3d = Ga('e_1 e_2 e_3', g=[1, 1, 1])
    assert o3d.g == Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


def test_constructor():
    """
    Test various multivector constructors
    """
    o3d, e_1, e_2, e_3 = Ga.build('e_1 e_2 e_3', g=[1,1,1])
    assert str(o3d.mv('a', 'scalar')) == 'a'
    assert str(o3d.mv('a', 'vector')) == 'a__1*e_1 + a__2*e_2 + a__3*e_3'
    assert str(o3d.mv('a', 'pseudo')) == 'a__123*e_1^e_2^e_3'
    assert str(o3d.mv('a', 'spinor')) == 'a + a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3'
    assert str(o3d.mv('a', 'mv')) == 'a + a__1*e_1 + a__2*e_2 + a__3*e_3 + a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3 + a__123*e_1^e_2^e_3'
    assert str(o3d.mv('a', 'bivector')) == 'a__12*e_1^e_2 + a__13*e_1^e_3 + a__23*e_2^e_3'


def test_basic_multivector_operations():

    g3d, ex, ey, ez = Ga.build('e*x|y|z')

    A = g3d.mv('A', 'mv')

    assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'

    X = g3d.mv('X', 'vector')
    Y = g3d.mv('Y', 'vector')

    assert str(X) == 'X__x*e_x + X__y*e_y + X__z*e_z'
    assert str(Y) == 'Y__x*e_x + Y__y*e_y + Y__z*e_z'

    assert str((X*Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z + (X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
    assert str((X ^ Y)) == '(X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
    assert str((X | Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z'

    g2d, ex, ey = Ga.build('e*x|y')

    X = g2d.mv('X', 'vector')
    A = g2d.mv('A', 'spinor')

    assert str(X) == 'X__x*e_x + X__y*e_y'
    assert str(A) == 'A + A__xy*e_x^e_y'

    assert str((X | A)) == '-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y)*e_x + A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y)*e_y'
    assert str((X < A)) == '(-(e_x.e_y)*A__xy*X__x - (e_y.e_y)*A__xy*X__y + A*X__x)*e_x + ((e_x.e_x)*A__xy*X__x + (e_x.e_y)*A__xy*X__y + A*X__y)*e_y'
    assert str((A > X)) == '((e_x.e_y)*A__xy*X__x + (e_y.e_y)*A__xy*X__y + A*X__x)*e_x + (-(e_x.e_x)*A__xy*X__x - (e_x.e_y)*A__xy*X__y + A*X__y)*e_y'

    o2d, ex, ey = Ga.build('e*x|y', g=[1, 1])

    X = o2d.mv('X', 'vector')
    A = o2d.mv('A', 'spinor')

    assert str(X) == 'X__x*e_x + X__y*e_y'
    assert str(A) == 'A + A__xy*e_x^e_y'

    assert str((X*A)) == '(A*X__x - A__xy*X__y)*e_x + (A*X__y + A__xy*X__x)*e_y'
    assert str((X | A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
    assert str((X < A)) == '(A*X__x - A__xy*X__y)*e_x + (A*X__y + A__xy*X__x)*e_y'
    assert str((X > A)) == 'A*X__x*e_x + A*X__y*e_y'

    assert str((A*X)) == '(A*X__x + A__xy*X__y)*e_x + (A*X__y - A__xy*X__x)*e_y'
    assert str((A | X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'
    assert str((A < X)) == 'A*X__x*e_x + A*X__y*e_y'
    assert str((A > X)) == '(A*X__x + A__xy*X__y)*e_x + (A*X__y - A__xy*X__x)*e_y'

    return

def test_check_generalized_BAC_CAB_formulas():

    g5d, a, b, c, d, e = Ga.build('a b c d e')

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

    X = (x, y, z) = symbols('x y z')
    o3d, ex, ey, ez = Ga.build('e_x e_y e_z', g=[1,1,1], coords=X)
    grad = o3d.grad

    f = o3d.mv('f', 'scalar', f=True)
    A = o3d.mv('A', 'vector', f=True)
    B = o3d.mv('B', 'bivector', f=True)
    C = o3d.mv('C', 'mv', f=True)

    assert str(f) == 'f'
    assert str(A) == 'A__x*e_x + A__y*e_y + A__z*e_z'
    assert str(B) == 'B__xy*e_x^e_y + B__xz*e_x^e_z + B__yz*e_y^e_z'
    assert str(C) == 'C + C__x*e_x + C__y*e_y + C__z*e_z + C__xy*e_x^e_y + C__xz*e_x^e_z + C__yz*e_y^e_z + C__xyz*e_x^e_y^e_z'

    assert str(grad*f) == 'D{x}f*e_x + D{y}f*e_y + D{z}f*e_z'
    assert str(grad | A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
    assert str(grad*A) == 'D{x}A__x + D{y}A__y + D{z}A__z + (-D{y}A__x + D{x}A__y)*e_x^e_y + (-D{z}A__x + D{x}A__z)*e_x^e_z + (-D{z}A__y + D{y}A__z)*e_y^e_z'

    assert str(-o3d.I()*(grad ^ A)) == '(-D{z}A__y + D{y}A__z)*e_x + (D{z}A__x - D{x}A__z)*e_y + (-D{y}A__x + D{x}A__y)*e_z'
    assert str(grad*B) == '(-D{y}B__xy - D{z}B__xz)*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z + (D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
    assert str(grad ^ B) == '(D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
    assert str(grad | B) == '(-D{y}B__xy - D{z}B__xz)*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'

    assert str(grad < A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
    assert str(grad > A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
    assert str(grad < B) == '(-D{y}B__xy - D{z}B__xz)*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'
    assert str(grad > B) == '0'
    assert str(grad < C) == 'D{x}C__x + D{y}C__y + D{z}C__z + (D{x}C - D{y}C__xy - D{z}C__xz)*e_x + (D{y}C + D{x}C__xy - D{z}C__yz)*e_y + (D{z}C + D{x}C__xz + D{y}C__yz)*e_z + D{z}C__xyz*e_x^e_y - D{y}C__xyz*e_x^e_z + D{x}C__xyz*e_y^e_z'
    assert str(grad > C) == 'D{x}C__x + D{y}C__y + D{z}C__z + D{x}C*e_x + D{y}C*e_y + D{z}C*e_z'

    return

def test_derivatives_in_spherical_coordinates():

    X = (r, th, phi) = symbols('r theta phi')
    sph3d, er, eth, ephi = Ga.build('e_r e_theta e_phi', g=[1,r**2,r**2*sin(th)**2], coords=X, norm=True)
    grad = sph3d.grad

    f = sph3d.mv('f', 'scalar', f=True)
    A = sph3d.mv('A', 'vector', f=True)
    B = sph3d.mv('B', 'bivector', f=True)

    assert str(f) == 'f'
    assert str(A) == 'A__r*e_r + A__theta*e_theta + A__phi*e_phi'
    assert str(B) == 'B__rtheta*e_r^e_theta + B__rphi*e_r^e_phi + B__thetaphi*e_theta^e_phi'

    assert str(grad*f) == 'D{r}f*e_r + D{theta}f*e_theta/r + D{phi}f*e_phi/(r*sin(theta))'
    assert str(grad | A) == 'D{r}A__r + (A__r + D{theta}A__theta)/r + (A__r*sin(theta) + A__theta*cos(theta) + D{phi}A__phi)/(r*sin(theta))'
    assert str(-sph3d.I()*(grad ^ A)) == '(A__phi/tan(theta) + D{theta}A__phi - D{phi}A__theta/sin(theta))*e_r/r + (-r*D{r}A__phi - A__phi + D{phi}A__r/sin(theta))*e_theta/r + (r*D{r}A__theta + A__theta - D{theta}A__r)*e_phi/r'
    assert str(grad ^ B) == '(r*D{r}B__thetaphi - B__rphi/tan(theta) + 2*B__thetaphi - D{theta}B__rphi + D{phi}B__rtheta/sin(theta))*e_r^e_theta^e_phi/r'

    return

def test_rounding_numerical_components():

    (o3d, ex, ey, ez) = Ga.build('e_x e_y e_z', g=[1,1,1])

    X = 1.2*ex + 2.34*ey + 0.555*ez
    Y = 0.333*ex + 4*ey + 5.3*ez

    assert str(X) == '1.2*e_x + 2.34*e_y + 0.555*e_z'
    assert str(Nga(X, 2)) == '1.2*e_x + 2.3*e_y + 0.55*e_z'
    assert str(X*Y) == '12.7011000000000 + 4.02078*e_x^e_y + 6.175185*e_x^e_z + 10.182*e_y^e_z'
    assert str(Nga(X*Y, 2)) == '13. + 4.0*e_x^e_y + 6.2*e_x^e_z + 10.0*e_y^e_z'

    return

def test_conformal_representations_of_circles_lines_spheres_and_planes():
    global n, nbar

    metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

    (cf3d, ex, ey, ez, n, nbar) = Ga.build('e_x e_y e_z n nbar', g=metric)

    x, y, z = symbols('x y z', real=True)

    e = n + nbar
    #conformal representation of points

    A = make_vector(ex, cf3d)
    B = make_vector(ey, cf3d)
    C = make_vector(-ex, cf3d)
    D = make_vector(ez, cf3d)
    X = make_vector(x*ex + y*ey +z*ez, cf3d)

    assert A == ex + (n - nbar)/S(2)
    assert B == ey + (n - nbar)/S(2)
    assert C == -ex + (n - nbar)/S(2)
    assert D == ez + (n - nbar)/S(2)
    assert X == x*ex + y*ey + z*ez + (x**2/2 + y**2/2 + z**2/2)*n - nbar/2

    assert A ^ B ^ C ^ X == -z*(ex^ey^ez^n) + z*(ex^ey^ez^nbar) + ((x**2 + y**2 + z**2 - S(1))/2)*(ex^ey^n^nbar)
    assert A ^ B ^ n ^ X == -z*(ex^ey^ez^n) + ((x + y - S(1))/2)*(ex^ey^n^nbar) + (z/2)*(ex^ez^n^nbar) - (z/2)*(ey^ez^n^nbar)
    assert A ^ B ^ C ^ D ^ X == ((-x**2 - y**2 - z**2 + S(1))/2)*ex^ey^ez^n^nbar
    assert A ^ B ^ n ^ D ^ X == ((-x - y - z + S(1))/2)*(ex^ey^ez^n^nbar)

    L = (A ^ B ^ e) ^ X

    assert L == -z*(ex^ey^ez^n) - z*(ex^ey^ez^nbar) + (-x**2/2 + x - y**2/2 + y - z**2/2 - S(1)/2)*(ex^ey^n^nbar) + z*(ex^ez^n^nbar) - z*(ey^ez^n^nbar)

    return

def test_properties_of_geometric_objects():

    metric = '# # # 0 0,' + \
             '# # # 0 0,' + \
             '# # # 0 0,' + \
             '0 0 0 0 2,' + \
             '0 0 0 2 0'

    (cf3d, p1, p2, p3, n, nbar) = Ga.build('p1 p2 p3 n nbar', g=metric)

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

    metric = '0 -1 #,' + \
             '-1 0 #,' + \
             '# # #'

    (cf1d, P1, P2, a) = Ga.build('P1 P2 a', g=metric)

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

def test_submanifolds():

    #Define spherical coordinate system in 3-d

    coords = (r, th, phi) = symbols('r,theta,phi', real=True)

    sp3d = Ga('e_r e_th e_ph', g=[1, r**2, r**2*sin(th)**2], coords=coords)
    (er, eth, ephi) = sp3d.mv()

    #Define coordinates for 2-d (u,v) and 1-d (s) manifolds

    u,v,s,alpha = symbols('u v s alpha',real=True)

    sub_coords = (u,v)

    smap = [1, u, v]  # Coordinate map for sphere of r = 1 in 3-d

    #Define unit sphere manifold

    sph2d = sp3d.sm(smap,sub_coords)
    (eu,ev) = sph2d.mv()

    #Define vector and vector field on unit sphere tangent space

    a = sph2d.mv('a','vector')
    b = sph2d.mv('b','vector')
    c = sph2d.mv('c','vector')
    f = sph2d.mv('f','vector',f=True)

    #Define directional derivative in direction a for unit sphere manifold

    dd = a|sph2d.grad

    assert str(dd) == 'a__u*D{u} + a__v*D{v}'
    assert str(dd * eu) == 'a__v*e_v/tan(u)'
    assert str(dd * ev) == '-a__v*sin(2*u)*e_u/2 + a__u*e_v/tan(u)'
    assert str(dd * f) == '(a__u*D{u}f__u - a__v*f__v*sin(2*u)/2 + a__v*D{v}f__u)*e_u + (a__u*f__v/tan(u) + a__u*D{u}f__v + a__v*f__u/tan(u) + a__v*D{v}f__v)*e_v'

    V = Mlt('V',sph2d,nargs=1,fct=True)
    T = Mlt('T',sph2d,nargs=2,fct=True)

    assert str(T.contract(1,2)) == 'a_1__u**2*D{u}^2T_uu + a_1__u**2*D{v}^2T_uu/sin(u)**2 + a_1__u*a_1__v*D{u}^2T_uv + a_1__u*a_1__v*D{u}^2T_vu + a_1__u*a_1__v*D{v}^2T_uv/sin(u)**2 + a_1__u*a_1__v*D{v}^2T_vu/sin(u)**2 + a_1__v**2*D{u}^2T_vv + a_1__v**2*D{v}^2T_vv/sin(u)**2'

    #Tensor Evaluation

    assert str(T(a,b)) == 'a__u*b__u*T_uu + a__u*b__v*T_uv + a__v*b__u*T_vu + a__v*b__v*T_vv'
    assert str(T(a,b+c).expand()) == 'a__u*b__u*T_uu + a__u*b__v*T_uv + a__u*c__u*T_uu + a__u*c__v*T_uv + a__v*b__u*T_vu + a__v*b__v*T_vv + a__v*c__u*T_vu + a__v*c__v*T_vv'
    assert str(T(a,alpha*b)) == 'a__u*alpha*b__u*T_uu + a__u*alpha*b__v*T_uv + a__v*alpha*b__u*T_vu + a__v*alpha*b__v*T_vv'

    #Geometric Derivative With Respect To Slot

    assert str(T.pdiff(1)) == '(a_1__u*a_2__u*D{u}T_uu + a_1__u*a_2__v*D{u}T_uv + a_1__v*a_2__u*D{u}T_vu + a_1__v*a_2__v*D{u}T_vv)*e_u + (a_1__u*a_2__u*D{v}T_uu + a_1__u*a_2__v*D{v}T_uv + a_1__v*a_2__u*D{v}T_vu + a_1__v*a_2__v*D{v}T_vv)*e_v/sin(u)**2'
    assert str(T.pdiff(2)) == '(a_1__u*a_2__u*D{u}T_uu + a_1__u*a_2__v*D{u}T_uv + a_1__v*a_2__u*D{u}T_vu + a_1__v*a_2__v*D{u}T_vv)*e_u + (a_1__u*a_2__u*D{v}T_uu + a_1__u*a_2__v*D{v}T_uv + a_1__v*a_2__u*D{v}T_vu + a_1__v*a_2__v*D{v}T_vv)*e_v/sin(u)**2'

    #Covariant Derivatives

    assert str(V.cderiv()) == 'a_2__u*(a_1__u*D{u}V_u + a_1__v*D{u}V_v) + a_2__v*(a_1__u*D{v}V_u + a_1__v*D{v}V_v)'

    DT = T.cderiv()

    assert str(DT) == 'a_3__u*(a_1__u*a_2__u*D{u}T_uu + a_1__u*a_2__v*D{u}T_uv + a_1__v*a_2__u*D{u}T_vu + a_1__v*a_2__v*D{u}T_vv) + a_3__v*(a_1__u*a_2__u*D{v}T_uu + a_1__u*a_2__v*D{v}T_uv + a_1__v*a_2__u*D{v}T_vu + a_1__v*a_2__v*D{v}T_vv)'

    #Define curve on unit sphere manifold

    us = Function('u__s')(s)
    vs = Function('v__s')(s)

    #Define 1-d submanifold on unit shpere manifold

    crv1d = sph2d.sm([us,vs],[s])

    (es,) = crv1d.mv()

    #1-D Manifold On Unit Sphere:

    assert str(crv1d.grad) == 'e_s*1/(sin(u__s)**2*D{s}v__s**2 + D{s}u__s**2)*D{s}'

    #Define scalar and vector fields on 1-d manifold tangent space

    g = crv1d.mv('g','scalar',f=True)
    h = crv1d.mv('h','vector',f=True)

    assert str(crv1d.grad * g) == 'D{s}g*e_s/(sin(u__s)**2*D{s}v__s**2 + D{s}u__s**2)'
    assert str(crv1d.grad | h) == 'D{s}h__s + h__s*sin(u__s)**2*D{s}v__s*D{s}^2v__s/(sin(u__s)**2*D{s}v__s**2 + D{s}u__s**2) + h__s*sin(2*u__s)*D{s}u__s*D{s}v__s**2/(2*(sin(u__s)**2*D{s}v__s**2 + D{s}u__s**2)) + h__s*D{s}u__s*D{s}^2u__s/(sin(u__s)**2*D{s}v__s**2 + D{s}u__s**2)'

    return

"""
test_differential_operators()
test_basic_multivector_operations()
test_check_generalized_BAC_CAB_formulas()
test_conformal_representations_of_circles_lines_spheres_and_planes()
test_constructor()
test_contraction()
test_derivative()
test_derivatives_in_rectangular_coordinates()
test_derivatives_in_spherical_coordinates()
test_extract_plane_and_line()
test_extracting_vectors_from_conformal_2_blade()
test_geometry()
test_metric()
test_metrics()
test_properties_of_geometric_objects()
test_reciprocal_frame()
test_rmul()
test_rounding_numerical_components()
test_substitution()
test_vector_extraction()
test_submanifolds()
"""
