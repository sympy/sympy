#!/usr/bin/python
#test_GA.py

"""
The reference D&L is "Geomertric Algebra for Physicists" by Doran and Lasenby
"""

try:
    import numpy
    disabled = False
except ImportError:
    #py.test will not execute any tests now
    disabled = True

if not disabled:
    from sympy.galgebra.GA import set_main, MV, make_symbols, types, ZERO, ONE, HALF
    import sympy
    from sympy import collect, sympify

    import sys
    set_main(sys.modules[__name__])

def F(x):
    """
    Conformal Mapping Function from 3D euclidian space to 5D conformal space
    where the images of all maps are null vectors.
    """
    Fx = HALF*((x*x)*n+2*x-nbar)
    #print 'F(x) =',Fx
    return(Fx)

def make_vector(a,n = 3):
    if type(a) == types.StringType:
        sym_str = ''
        for i in range(n):
            sym_str += a+str(i)+' '
        sym_lst = make_symbols(sym_str)
        sym_lst.append(ZERO)
        sym_lst.append(ZERO)
        a = MV(sym_lst,'vector')
    return(F(a))

def test_rmul():
    """
    Test for communitive scalar multiplication.  Leftover from when sympy and
    numpy were not working together and __mul__ and __rmul__ would not give the
    same answer.
    """
    MV.setup('x y z')
    make_symbols('a b c')
    assert 5*x == x*5
    assert HALF*x == x*HALF
    assert a*x == x*a

def test_noneuclidian():
    """
    Test of complex geometric algebra manipulation to derive distance function
    for 2-D hyperbolic non-euclidian space.  See D&L Section 10.6.2
    """
    metric = '0 # #,'+ \
             '# 0 #,'+ \
             '# # 1,'
    MV.setup('X Y e',metric,debug=0)
    MV.set_str_format(1)
    L = X^Y^e
    B = L*e
    Bsq = (B*B)()
    BeBr =B*e*B.rev()
    (s,c,Binv,M,S,C,alpha) = sympy.symbols('s','c','Binv','M','S','C','alpha')
    Bhat = Binv*B # Normalize translation generator
    R = c+s*Bhat # Rotor R = exp(alpha*Bhat/2)
    Z = R*X*R.rev()
    Z.expand()
    Z.collect([Binv,s,c,XdotY])
    W = Z|Y
    W.expand()
    W.collect([s*Binv])
    M = 1/Bsq
    W.subs(Binv**2,M)
    W.simplify()
    Bmag = sympy.sqrt(XdotY**2-2*XdotY*Xdote*Ydote)
    W.collect([Binv*c*s,XdotY])
    W.subs(2*XdotY**2-4*XdotY*Xdote*Ydote,2/(Binv**2))
    W.subs(2*c*s,S)
    W.subs(c**2,(C+1)/2)
    W.subs(s**2,(C-1)/2)
    W.simplify()
    W.subs(1/Binv,Bmag)

    W = W().expand()
    #print '(R*X*R.rev()).Y =',W
    Wd = collect(W,[C,S],exact=True,evaluate=False)
    #print 'Wd =',Wd
    Wd_1 = Wd[ONE]
    Wd_C = Wd[C]
    Wd_S = Wd[S]
    #print '|B| =',Bmag
    Wd_1 = Wd_1.subs(Bmag,1/Binv)
    Wd_C = Wd_C.subs(Bmag,1/Binv)
    Wd_S = Wd_S.subs(Bmag,1/Binv)
    #print 'Wd[ONE] =',Wd_1
    #print 'Wd[C] =',Wd_C
    #print 'Wd[S] =',Wd_S
    lhs = Wd_1+Wd_C*C
    rhs = -Wd_S*S
    lhs = lhs**2
    rhs = rhs**2
    W = (lhs-rhs).expand()
    W = (W.subs(1/Binv**2,Bmag**2)).expand()
    #print 'W =',W
    W = (W.subs(S**2,C**2-1)).expand()
    W = collect(W,[C**2,C],evaluate=False)
    #print 'W =',W
    a = W[C**2]
    b = W[C]
    c = W[ONE]
    #print 'a =',a
    #print 'b =',b
    #print 'c =',c
    D = (b**2-4*a*c).expand()
    #print 'Setting to 0 and solving for C gives:'
    #print 'Descriminant D = b^2-4*a*c =',D
    C = (-b/(2*a)).expand()
    #print 'C = cosh(alpha) = -b/(2*a) =',C

    #cosh(alpha) = 1-X.Y/((X.e)(Y.e))
    #alpha is noneuclidian distance
    assert C == 1-XdotY/(Xdote*Ydote)

def test_reciprocal_frame():
    """
    Test of fromula for general reciprocal frame of three vectors.
    Let three independent vectors be e1, e2, and e3. The reciprocal
    vectors E1, E2, and E3 obey the relations:

    e_i.E_j = delta_ij*(e1^e2^e3)**2
    """
    metric = '1 # #,'+ \
             '# 1 #,'+ \
             '# # 1,'

    MV.setup('e1 e2 e3',metric)
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

def test_vector_extraction():
    """
    Show that conformal bivector encodes two points. See D&L Section 10.4.1
    """
    metric = ' 0 -1 #,'+ \
             '-1  0 #,'+ \
             ' #  # #,'

    MV.setup('P1 P2 a',metric)
    """
    P1 and P2 are null vectors and hence encode points in conformal space.
    Show that P1 and P2 can be extracted from the bivector B = P1^P2. a is a
    third vector in the conformal space with a.B not 0.
    """
    ZERO_MV = MV()
    B = P1^P2
    Bsq = B*B
    ap = a-(a^B)*B
    Ap = ap+ap*B
    Am = ap-ap*B
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
    assert Ap2 == ZERO_MV
    assert Am2 == ZERO_MV

def test_geometry():
    """
    Test conformal geometric description of circles, lines, spheres, and planes.
    """
    metric = '1 0 0 0 0,'+ \
             '0 1 0 0 0,'+ \
             '0 0 1 0 0,'+ \
             '0 0 0 0 2,'+ \
             '0 0 0 2 0'

    MV.setup('e0 e1 e2 n nbar',metric,debug=0)
    e = n+nbar
    #conformal representation of points
    ZERO_MV = MV()

    A = make_vector(e0)    # point a = (1,0,0)  A = F(a)
    B = make_vector(e1)    # point b = (0,1,0)  B = F(b)
    C = make_vector(-1*e0) # point c = (-1,0,0) C = F(c)
    D = make_vector(e2)    # point d = (0,0,1)  D = F(d)
    X = make_vector('x',3)

    Circle = A^B^C^X
    Line = A^B^n^X
    Sphere = A^B^C^D^X
    Plane = A^B^n^D^X

    #Circle through a, b, and c
    Circle_test = -x2*(e0^e1^e2^n)+x2*(e0^e1^e2^nbar)+HALF*(-1+x0**2+x1**2+x2**2)*(e0^e1^n^nbar)
    diff = Circle-Circle_test
    diff.compact()
    assert diff == ZERO_MV

    #Line through a and b
    Line_test = -x2*(e0^e1^e2^n)+HALF*(-1+x0+x1)*(e0^e1^n^nbar)+(HALF*x2)*(e0^e2^n^nbar)+\
                (-HALF*x2)*(e1^e2^n^nbar)
    diff = Line-Line_test
    diff.compact()
    assert diff == ZERO_MV

    #Sphere through a, b, c, and d
    Sphere_test = HALF*(1-x0**2-x1**2-x2**2)*(e0^e1^e2^n^nbar)
    diff = Sphere-Sphere_test
    diff.compact()
    assert diff == ZERO_MV

    #Plane through a, b, and d
    Plane_test = HALF*(1-x0-x1-x2)*(e0^e1^e2^n^nbar)
    diff = Plane-Plane_test
    diff.compact()
    assert diff == ZERO_MV

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

    MV.setup('p1 p2 p3 n nbar',metric,debug=0)
    MV.set_str_format(1)

    ZERO_MV = MV()

    P1 = F(p1)
    P2 = F(p2)
    P3 = F(p3)

    #Line through p1 and p2
    L = P1^P2^n
    delta = (L|n)|nbar
    delta_test = 2*p1-2*p2
    diff = delta-delta_test
    diff.compact()
    assert diff == ZERO_MV

    #Plane through p1, p2, and p3
    C = P1^P2^P3
    delta = ((C^n)|n)|nbar
    delta_test = 2*(p1^p2)-2*(p1^p3)+2*(p2^p3)
    diff = delta-delta_test
    diff.compact()
    assert diff == ZERO_MV
