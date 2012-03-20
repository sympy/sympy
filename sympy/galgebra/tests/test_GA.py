#!/usr/bin/python
#test_GA.py

"""
The reference D&L is "Geometric Algebra for Physicists" by Doran and Lasenby
"""

import sys

from sympy.external import import_module
numpy = import_module('numpy')
if not numpy:
    disabled = True
else:
    sys.path.append('../')
    from sympy.galgebra.GA import MV, ZERO, HALF, S
    import sympy

def F(x, n, nbar):
    """
    Conformal Mapping Function from 3D Euclidean space to 5D conformal space
    where the images of all maps are null vectors.
    """
    Fx = HALF*((x*x)*n+2*x-nbar)
    #print 'F(x) =',Fx
    return(Fx)

def test_rmul():
    """
    Test for commutative scalar multiplication.  Leftover from when sympy and
    numpy were not working together and __mul__ and __rmul__ would not give the
    same answer.
    """
    x,y,z = MV.setup('x y z')
    a,b,c = sympy.symbols('a b c')
    assert 5*x == x*5
    assert HALF*x == x*HALF
    assert a*x == x*a

def test_contraction():
    """
    Test for inner product and left and right contraction
    """

    e_1,e_2,e_3 = MV.setup('e_1 e_2 e_3','1 0 0, 0 1 0, 0 0 1',offset=1)

    assert ((e_1^e_3)|e_1) == -e_3
    assert ((e_1^e_3)>e_1) == -e_3
    assert (e_1|(e_1^e_3)) == e_3
    assert (e_1<(e_1^e_3)) == e_3
    assert ((e_1^e_3)<e_1) == 0
    assert (e_1>(e_1^e_3)) == 0

def test_substitution():

    e_x,e_y,e_z = MV.setup('e_x e_y e_z','1 0 0, 0 1 0, 0 0 1',offset=1)
    x,y,z = sympy.symbols('x y z')

    X = x*e_x+y*e_y+z*e_z
    Y = X.subs([(x,2),(y,3),(z,4)])
    assert Y == 2*e_x+3*e_y+4*e_z


def test_vector_extraction():
    """
    Show that conformal bivector encodes two points. See D&L Section 10.4.1
    """
    metric = ' 0 -1 #,'+ \
             '-1  0 #,'+ \
             ' #  # #,'

    P1,P2,a = MV.setup('P1 P2 a',metric)
    """
    P1 and P2 are null vectors and hence encode points in conformal space.
    Show that P1 and P2 can be extracted from the bivector B = P1^P2. a is a
    third vector in the conformal space with a.B not 0.
    """
    B = P1^P2
    Bsq = B*B
    ap = a-(a^B)*B
    Ap = ap+ap*B
    Am = ap-ap*B
    P1dota = sympy.Symbol('(P1.a)')
    P2dota = sympy.Symbol('(P2.a)')
    Ap_test = (-2*P2dota)*P1
    Am_test = (-2*P1dota)*P2
    Ap.compact()
    Am.compact()
    Ap_test.compact()
    Am_test.compact()
    assert Ap == Ap_test
    assert Am == Am_test
    Ap2 = Ap*Ap
    Am2 = Am*Am
    Ap2.compact()
    Am2.compact()
    assert Ap2 == ZERO
    assert Am2 == ZERO

def test_geometry():
    """
    Test conformal geometric description of circles, lines, spheres, and planes.
    """
    metric = '1 0 0 0 0,'+ \
             '0 1 0 0 0,'+ \
             '0 0 1 0 0,'+ \
             '0 0 0 0 2,'+ \
             '0 0 0 2 0'

    e0,e1,e2,n,nbar = MV.setup('e0 e1 e2 n nbar',metric,debug=0)
    e = n+nbar
    #conformal representation of points

    A = F(e0,n,nbar)    # point a = (1,0,0)  A = F(a)
    B = F(e1,n,nbar)    # point b = (0,1,0)  B = F(b)
    C = F(-1*e0,n,nbar) # point c = (-1,0,0) C = F(c)
    D = F(e2,n,nbar)    # point d = (0,0,1)  D = F(d)
    x0,x1,x2 = sympy.symbols('x0 x1 x2')
    X = F(MV([x0,x1,x2],'vector'),n,nbar)

    Circle = A^B^C^X
    Line = A^B^n^X
    Sphere = A^B^C^D^X
    Plane = A^B^n^D^X

    #Circle through a, b, and c
    Circle_test = -x2*(e0^e1^e2^n)+x2*(e0^e1^e2^nbar)+HALF*(-1+x0**2+x1**2+x2**2)*(e0^e1^n^nbar)
    diff = Circle-Circle_test
    diff.compact()
    assert diff == ZERO

    #Line through a and b
    Line_test = -x2*(e0^e1^e2^n)+HALF*(-1+x0+x1)*(e0^e1^n^nbar)+(HALF*x2)*(e0^e2^n^nbar)+\
                (-HALF*x2)*(e1^e2^n^nbar)
    diff = Line-Line_test
    diff.compact()
    assert diff == ZERO

    #Sphere through a, b, c, and d
    Sphere_test = HALF*(1-x0**2-x1**2-x2**2)*(e0^e1^e2^n^nbar)
    diff = Sphere-Sphere_test
    diff.compact()
    assert diff == ZERO

    #Plane through a, b, and d
    Plane_test = HALF*(1-x0-x1-x2)*(e0^e1^e2^n^nbar)
    diff = Plane-Plane_test
    diff.compact()
    assert diff == ZERO

def test_extract_plane_and_line():
    """
    Show that conformal trivector encodes planes and lines. See D&L section
    10.4.2
    """
    metric = '# # # 0 0,'+ \
             '# # # 0 0,'+ \
             '# # # 0 0,'+ \
             '0 0 0 0 2,'+ \
             '0 0 0 2 0'

    p1,p2,p3,n,nbar = MV.setup('p1 p2 p3 n nbar',metric,debug=0)
    MV.set_str_format(1)

    P1 = F(p1,n,nbar)
    P2 = F(p2,n,nbar)
    P3 = F(p3,n,nbar)

    #Line through p1 and p2
    L = P1^P2^n
    delta = (L|n)|nbar
    delta_test = 2*p1-2*p2
    diff = delta-delta_test
    diff.compact()
    assert diff == ZERO

    #Plane through p1, p2, and p3
    C = P1^P2^P3
    delta = ((C^n)|n)|nbar
    delta_test = 2*(p1^p2)-2*(p1^p3)+2*(p2^p3)
    diff = delta-delta_test
    diff.compact()
    assert diff == ZERO

def test_reciprocal_frame():
    """
    Test of formula for general reciprocal frame of three vectors.
    Let three independent vectors be e1, e2, and e3. The reciprocal
    vectors E1, E2, and E3 obey the relations:

    e_i.E_j = delta_ij*(e1^e2^e3)**2
    """
    metric = '1 # #,'+ \
             '# 1 #,'+ \
             '# # 1,'

    e1,e2,e3 = MV.setup('e1 e2 e3',metric)
    E = e1^e2^e3
    Esq = (E*E)()
    Esq_inv = 1/Esq
    E1 = (e2^e3)*E
    E2 = (-1)*(e1^e3)*E
    E3 = (e1^e2)*E
    w = (E1|e2)
    w.collect(MV.g)
    w = w().expand()
    w = (E1|e3)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E2|e1)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E2|e3)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E3|e1)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E3|e2)
    w.collect(MV.g)
    w = w().expand()
    assert w == 0
    w = (E1|e1)
    w = w().expand()
    Esq = Esq.expand()
    assert w/Esq == 1
    w = (E2|e2)
    w = w().expand()
    assert w/Esq == 1
    w = (E3|e3)
    w = w().expand()
    assert w/Esq == 1

def test_derivative():
    coords = x,y,z = sympy.symbols('x y z')
    e_x,e_y,e_z = MV.setup('e','1 0 0, 0 1 0, 0 0 1',coords=coords)
    X = x*e_x+y*e_y+z*e_z
    a = MV('a','vector')

    assert ((X|a).grad()) == a
    assert ((X*X).grad()) == 2*X
    assert (X*X*X).grad() == 5*X*X
    assert X.grad_int() == 3

def test_str():
    e_1,e_2,e_3 = MV.setup('e_1 e_2 e_3','1 0 0, 0 1 0, 0 0 1')

    X = MV('x')
    assert str(X) == 'x+x__0*e_1+x__1*e_2+x__2*e_3+x__01*e_1e_2+x__02*e_1e_3+x__12*e_2e_3+x__012*e_1e_2e_3'
    Y = MV('y','spinor')
    assert str(Y) == 'y+y__01*e_1e_2+y__02*e_1e_3+y__12*e_2e_3'
    Z = X+Y
    assert str(Z) == 'x+y+x__0*e_1+x__1*e_2+x__2*e_3+(x__01+y__01)*e_1e_2+(x__02+y__02)*e_1e_3+(x__12+y__12)*e_2e_3+x__012*e_1e_2e_3'
    assert str(e_1|e_1) == '1'

def test_metric():
    MV.setup('e_1 e_2 e_3','[1,1,1]')
    assert str(MV.metric) == '[[1 0 0]\n [0 1 0]\n [0 0 1]]'

def test_constructor():
    """
    Test various multivector constructors
    """
    e_1,e_2,e_3 = MV.setup('e_1 e_2 e_3','[1,1,1]')
    x = sympy.symbols('x')
    assert str(S(1)) == '1'
    assert str(S(x)) == 'x'
    assert str(MV('a','scalar')) == 'a'
    assert str(MV('a','vector')) == 'a__0*e_1+a__1*e_2+a__2*e_3'
    assert str(MV('a','pseudo')) == 'a*e_1e_2e_3'
    assert str(MV('a','spinor')) == 'a+a__01*e_1e_2+a__02*e_1e_3+a__12*e_2e_3'
    assert str(MV('a')) == 'a+a__0*e_1+a__1*e_2+a__2*e_3+a__01*e_1e_2+a__02*e_1e_3+a__12*e_2e_3+a__012*e_1e_2e_3'
    assert str(MV([2,'a'],'grade')) == 'a__01*e_1e_2+a__02*e_1e_3+a__12*e_2e_3'
    assert str(MV('a','grade2')) == 'a__01*e_1e_2+a__02*e_1e_3+a__12*e_2e_3'

def test__print_Mul_Add():
    from sympy.galgebra.latex_ex import LatexPrinter
    from sympy import symbols
    n, m = symbols('n,m', negative=True)
    z = symbols('z')
    l = LatexPrinter()
    assert l._print_Mul(n*m) == 'm n'
    assert l._print_Mul(-2*m) == '- 2 m'
    assert l._print_Mul(2*m) == '2 m'
    assert l._print_Add(-5 + 4*z) == '-5 + 4 z'
    assert l._print_Add(-5 - 4*z) == '-5 - 4 z'
    assert l._print_Add(n - 2) == '-2 + n'
