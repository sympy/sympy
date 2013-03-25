#!/usr/bin/python
import sys

from sympy import symbols,sin,cos
from sympy.GA.GA import *
from sympy.GA.GAPrint import *


def test_basic_multivector_operations():
    (ex,ey,ez) = MV.setup('e*x|y|z')

    A = MV('A','mv')

    assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'
    assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'
    assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'

    X = MV('X','vector')
    Y = MV('Y','vector')


    assert str(X) == 'X__x*e_x + X__y*e_y + X__z*e_z'
    assert str(Y) == 'Y__x*e_x + Y__y*e_y + Y__z*e_z'

    assert str((X*Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z + (X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
    assert str((X^Y)) == '(X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
    assert str((X|Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z'

    (ex,ey) = MV.setup('e*x|y')


    X = MV('X','vector')
    A = MV('A','spinor')

    assert str(X) == 'X__x*e_x + X__y*e_y'
    assert str(A) == 'A + A__xy*e_x^e_y'

    assert str((X|A)) == '(-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y))*e_x + (A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y))*e_y'
    assert str((X<A)) == '(-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y))*e_x + (A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y))*e_y'
    assert str((A>X)) == '(A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y))*e_x + (-A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y))*e_y'

    (ex,ey) = MV.setup('e*x|y',metric='[1,1]')


    X = MV('X','vector')
    A = MV('A','spinor')

    assert str(X) == 'X__x*e_x + X__y*e_y'
    assert str(A) == 'A + A__xy*e_x^e_y'

    assert str((X*A)) == '(A*X__x - A__xy*X__y)*e_x + (A*X__y + A__xy*X__x)*e_y'
    assert str((X|A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
    assert str((X<A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
    assert str((X>A)) == 'A*X__x*e_x + A*X__y*e_y'

    assert str((A*X)) == '(A*X__x + A__xy*X__y)*e_x + (A*X__y - A__xy*X__x)*e_y'
    assert str((A|X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'
    assert str((A<X)) == 'A*X__x*e_x + A*X__y*e_y'
    assert str((A>X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'
    return

def test_check_generalized_BAC_CAB_formulas():

    (a,b,c,d,e) = MV.setup('a b c d e')


    assert str(a|(b*c)) == '-(a.c)*b + (a.b)*c'
    assert str(a|(b^c)) == '-(a.c)*b + (a.b)*c'
    assert str(a|(b^c^d)) == '(a.d)*b^c - (a.c)*b^d + (a.b)*c^d'
    assert str((a|(b^c))+(c|(a^b))+(b|(c^a))) == '0'
    assert str(a*(b^c)-b*(a^c)+c*(a^b)) == '3*a^b^c'
    assert str(a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)) == '4*a^b^c^d'
    assert str((a^b)|(c^d)) == '-(a.c)*(b.d) + (a.d)*(b.c)'
    assert str(((a^b)|c)|d) == '-(a.c)*(b.d) + (a.d)*(b.c)'
    assert str(Com(a^b,c^d)) == '-(b.d)*a^c + (b.c)*a^d + (a.d)*b^c - (a.c)*b^d'
    assert str((a|(b^c))|(d^e)) == '(-(a.b)*(c.e) + (a.c)*(b.e))*d + ((a.b)*(c.d) - (a.c)*(b.d))*e'

    return

def test_derivatives_in_rectangular_coordinates():

    X = (x,y,z) = symbols('x y z')
    (ex,ey,ez,grad) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=X)

    f = MV('f','scalar',fct=True)
    A = MV('A','vector',fct=True)
    B = MV('B','grade2',fct=True)
    C = MV('C','mv',fct=True)
    assert str(f) == 'f'
    assert str(A) == 'A__x*e_x + A__y*e_y + A__z*e_z'
    assert str(B) == 'B__xy*e_x^e_y + B__xz*e_x^e_z + B__yz*e_y^e_z'
    assert str(C) == 'C + C__x*e_x + C__y*e_y + C__z*e_z + C__xy*e_x^e_y + C__xz*e_x^e_z + C__yz*e_y^e_z + C__xyz*e_x^e_y^e_z'

    assert str(grad*f) == 'D{x}f*e_x + D{y}f*e_y + D{z}f*e_z'
    assert str(grad|A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
    assert str(grad*A) == 'D{x}A__x + D{y}A__y + D{z}A__z + (-D{y}A__x + D{x}A__y)*e_x^e_y + (-D{z}A__x + D{x}A__z)*e_x^e_z + (-D{z}A__y + D{y}A__z)*e_y^e_z'

    assert str(-MV.I*(grad^A)) == '(-D{z}A__y + D{y}A__z)*e_x + (D{z}A__x - D{x}A__z)*e_y + (-D{y}A__x + D{x}A__y)*e_z'
    assert str(grad*B) == '(-(D{y}B__xy + D{z}B__xz))*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z + (D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
    assert str(grad^B) == '(D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
    assert str(grad|B) == '(-(D{y}B__xy + D{z}B__xz))*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'

    assert str(grad<A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
    assert str(grad>A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
    assert str(grad<B) == '(-(D{y}B__xy + D{z}B__xz))*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'
    assert str(grad>B) == '0'
    assert str(grad<C) == 'D{x}C__x + D{y}C__y + D{z}C__z + (-(D{y}C__xy + D{z}C__xz))*e_x + (D{x}C__xy - D{z}C__yz)*e_y + (D{x}C__xz + D{y}C__yz)*e_z + D{z}C__xyz*e_x^e_y - D{y}C__xyz*e_x^e_z + D{x}C__xyz*e_y^e_z'
    assert str(grad>C) == 'D{x}C__x + D{y}C__y + D{z}C__z + D{x}C*e_x + D{y}C*e_y + D{z}C*e_z'

    return

def test_derivatives_in_spherical_coordinates():

    X = (r,th,phi) = symbols('r theta phi')
    curv = [[r*cos(phi)*sin(th),r*sin(phi)*sin(th),r*cos(th)],[1,r,r*sin(th)]]
    (er,eth,ephi,grad) = MV.setup('e_r e_theta e_phi',metric='[1,1,1]',coords=X,curv=curv)

    f = MV('f','scalar',fct=True)
    A = MV('A','vector',fct=True)
    B = MV('B','grade2',fct=True)

    assert str(f) == 'f'
    assert str(A) == 'A__r*e_r + A__theta*e_theta + A__phi*e_phi'
    assert str(B) == 'B__rtheta*e_r^e_theta + B__rphi*e_r^e_phi + B__thetaphi*e_theta^e_phi'

    assert str(grad*f) == 'D{r}f*e_r + D{theta}f/r*e_theta + D{phi}f/(r*sin(theta))*e_phi'
    assert str(grad|A) == 'D{r}A__r + 2*A__r/r + A__theta*cos(theta)/(r*sin(theta)) + D{theta}A__theta/r + D{phi}A__phi/(r*sin(theta))'
    assert str(-MV.I*(grad^A)) == '((A__phi*cos(theta)/sin(theta) + D{theta}A__phi - D{phi}A__theta/sin(theta))/r)*e_r + (-D{r}A__phi - A__phi/r + D{phi}A__r/(r*sin(theta)))*e_theta + (D{r}A__theta + A__theta/r - D{theta}A__r/r)*e_phi'
    assert str(grad^B) == '(D{r}B__thetaphi - B__rphi*cos(theta)/(r*sin(theta)) + 2*B__thetaphi/r - D{theta}B__rphi/r + D{phi}B__rtheta/(r*sin(theta)))*e_r^e_theta^e_phi'
    return

def test_rounding_numerical_components():

    (ex,ey,ez) = MV.setup('e_x e_y e_z',metric='[1,1,1]')

    X = 1.2*ex+2.34*ey+0.555*ez
    Y = 0.333*ex+4*ey+5.3*ez

    assert str(X) == '1.20000000000000*e_x + 2.34000000000000*e_y + 0.555000000000000*e_z'
    assert str(Nga(X,2)) == '1.2*e_x + 2.3*e_y + 0.55*e_z'
    assert str(X*Y) == '12.7011000000000 + 4.02078000000000*e_x^e_y + 6.17518500000000*e_x^e_z + 10.1820000000000*e_y^e_z'
    assert str(Nga(X*Y,2)) == '13. + 4.0*e_x^e_y + 6.2*e_x^e_z + 10.*e_y^e_z'
    return

def test_noneuclidian_distance_calculation():
    from sympy import solve,sqrt

    metric = '0 # #,# 0 #,# # 1'
    (X,Y,e) = MV.setup('X Y e',metric)


    assert str((X^Y)*(X^Y)) == '(X.Y)**2'

    L = X^Y^e
    B = L*e
    assert str(B) == 'X^Y - (Y.e)*X^e + (X.e)*Y^e'
    Bsq = B*B
    assert str(Bsq) == '(X.Y)*((X.Y) - 2*(X.e)*(Y.e))'
    Bsq = Bsq.scalar()
    assert str(B) == 'X^Y - (Y.e)*X^e + (X.e)*Y^e'

    BeBr =B*e*B.rev()
    assert str(BeBr) == '((X.Y)*(-(X.Y) + 2*(X.e)*(Y.e)))*e'
    assert str(B*B) == '(X.Y)*((X.Y) - 2*(X.e)*(Y.e))'
    assert str(L*L) == '(X.Y)*((X.Y) - 2*(X.e)*(Y.e))'
    (s,c,Binv,M,S,C,alpha,XdotY,Xdote,Ydote) = symbols('s c (1/B) M S C alpha (X.Y) (X.e) (Y.e)')

    Bhat = Binv*B
    R = c+s*Bhat
    assert str(R) == 'c + (1/B)*s*X^Y - (1/B)*(Y.e)*s*X^e + (1/B)*(X.e)*s*Y^e'

    Z = R*X*R.rev()
    Z.obj = expand(Z.obj)
    Z.obj = Z.obj.collect([Binv,s,c,XdotY])
    assert str(Z) == '((1/B)**2*(X.Y)**2*s**2 - 2*(1/B)**2*(X.Y)*(X.e)*(Y.e)*s**2 + 2*(1/B)*(X.Y)*c*s - 2*(1/B)*(X.e)*(Y.e)*c*s + c**2)*X + 2*(1/B)*(X.e)**2*c*s*Y + (2*(1/B)*(X.Y)*(X.e)*s*(-(1/B)*(X.Y)*s + 2*(1/B)*(X.e)*(Y.e)*s - c))*e'
    W = Z|Y
    # From this point forward all calculations are with sympy scalars
    W = W.scalar()
    assert str(W) == '(1/B)**2*(X.Y)**3*s**2 - 4*(1/B)**2*(X.Y)**2*(X.e)*(Y.e)*s**2 + 4*(1/B)**2*(X.Y)*(X.e)**2*(Y.e)**2*s**2 + 2*(1/B)*(X.Y)**2*c*s - 4*(1/B)*(X.Y)*(X.e)*(Y.e)*c*s + (X.Y)*c**2'
    W = expand(W)
    W = simplify(W)
    W = W.collect([s*Binv])

    M = 1/Bsq
    W = W.subs(Binv**2,M)
    W = simplify(W)
    Bmag = sqrt(XdotY**2-2*XdotY*Xdote*Ydote)
    W = W.collect([Binv*c*s,XdotY])

    #Double angle substitutions

    W = W.subs(2*XdotY**2-4*XdotY*Xdote*Ydote,2/(Binv**2))
    W = W.subs(2*c*s,S)
    W = W.subs(c**2,(C+1)/2)
    W = W.subs(s**2,(C-1)/2)
    W = simplify(W)
    W = W.subs(1/Binv,Bmag)
    W = expand(W)


    assert str(W) == '(X.Y)*C - (X.e)*(Y.e)*C + (X.e)*(Y.e) + S*sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'

    Wd = collect(W,[C,S],exact=True,evaluate=False)

    Wd_1 = Wd[ONE]
    Wd_C = Wd[C]
    Wd_S = Wd[S]

    assert str(Wd_1) == '(X.e)*(Y.e)'
    assert str(Wd_C) == '(X.Y) - (X.e)*(Y.e)'
    assert str(Wd_S) == 'sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'

    assert str(Bmag) == 'sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'
    Wd_1 = Wd_1.subs(Bmag,1/Binv)
    Wd_C = Wd_C.subs(Bmag,1/Binv)
    Wd_S = Wd_S.subs(Bmag,1/Binv)

    lhs = Wd_1+Wd_C*C
    rhs = -Wd_S*S
    lhs = lhs**2
    rhs = rhs**2
    W = expand(lhs-rhs)
    W = expand(W.subs(1/Binv**2,Bmag**2))
    W = expand(W.subs(S**2,C**2-1))
    W = W.collect([C,C**2],evaluate=False)

    a = simplify(W[C**2])
    b = simplify(W[C])
    c = simplify(W[ONE])


    assert str(a) == '(X.e)**2*(Y.e)**2'
    assert str(b) == '2*(X.e)*(Y.e)*((X.Y) - (X.e)*(Y.e))'
    assert str(c) == '(X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e) + (X.e)**2*(Y.e)**2'

    x = Symbol('x')
    C =  solve(a*x**2+b*x+c,x)[0]
    assert str(expand(simplify(expand(C)))) == '-(X.Y)/((X.e)*(Y.e)) + 1'
    return

HALF = Rational(1,2)

def F(x):
    global n,nbar
    Fx = HALF*((x*x)*n+2*x-nbar)
    return(Fx)

def make_vector(a,n = 3):
    if isinstance(a,str):
        sym_str = ''
        for i in range(n):
            sym_str += a+str(i+1)+' '
        sym_lst = list(symbols(sym_str))
        sym_lst.append(ZERO)
        sym_lst.append(ZERO)
        a = MV(sym_lst,'vector')
    return(F(a))

def test_conformal_representations_of_circles_lines_spheres_and_planes():
    global n,nbar

    metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

    (e1,e2,e3,n,nbar) = MV.setup('e_1 e_2 e_3 n nbar',metric)


    e = n+nbar
    #conformal representation of points

    A = make_vector(e1)
    B = make_vector(e2)
    C = make_vector(-e1)
    D = make_vector(e3)
    X = make_vector('x',3)

    assert str(A) == 'e_1 + 1/2*n - 1/2*nbar'
    assert str(B) == 'e_2 + 1/2*n - 1/2*nbar'
    assert str(C) == '-e_1 + 1/2*n - 1/2*nbar'
    assert str(D) == 'e_3 + 1/2*n - 1/2*nbar'
    assert str(X) == 'x1*e_1 + x2*e_2 + x3*e_3 + ((x1**2 + x2**2 + x3**2)/2)*n - 1/2*nbar'

    assert str((A^B^C^X)) == '-x3*e_1^e_2^e_3^n + x3*e_1^e_2^e_3^nbar + ((x1**2 + x2**2 + x3**2 - 1)/2)*e_1^e_2^n^nbar'
    assert str((A^B^n^X)) == '-x3*e_1^e_2^e_3^n + ((x1 + x2 - 1)/2)*e_1^e_2^n^nbar + x3/2*e_1^e_3^n^nbar - x3/2*e_2^e_3^n^nbar'
    assert str((((A^B)^C)^D)^X) == '((-x1**2 - x2**2 - x3**2 + 1)/2)*e_1^e_2^e_3^n^nbar'
    assert str((A^B^n^D^X)) == '((-x1 - x2 - x3 + 1)/2)*e_1^e_2^e_3^n^nbar'

    L = (A^B^e)^X

    assert str(L) == '-x3*e_1^e_2^e_3^n - x3*e_1^e_2^e_3^nbar + (-x1**2/2 + x1 - x2**2/2 + x2 - x3**2/2 - 1/2)*e_1^e_2^n^nbar + x3*e_1^e_3^n^nbar - x3*e_2^e_3^n^nbar'
    return

def test_properties_of_geometric_objects():

    metric = '# # # 0 0,'+ \
             '# # # 0 0,'+ \
             '# # # 0 0,'+ \
             '0 0 0 0 2,'+ \
             '0 0 0 2 0'

    (p1,p2,p3,n,nbar) = MV.setup('p1 p2 p3 n nbar',metric)


    P1 = F(p1)
    P2 = F(p2)
    P3 = F(p3)


    L = P1^P2^n
    delta = (L|n)|nbar
    assert str(delta) == '2*p1 - 2*p2'


    C = P1^P2^P3
    delta = ((C^n)|n)|nbar
    assert str(delta) == '2*p1^p2 - 2*p1^p3 + 2*p2^p3'
    assert str((p2-p1)^(p3-p1)) == 'p1^p2 - p1^p3 + p2^p3'

def test_extracting_vectors_from_conformal_2_blade():

    metric = ' 0 -1 #,'+ \
             '-1  0 #,'+ \
             ' #  # #,'

    (P1,P2,a) = MV.setup('P1 P2 a',metric)


    B = P1^P2
    Bsq = B*B
    assert str(Bsq) == '1'
    ap = a-(a^B)*B
    assert str(ap) == '-(P2.a)*P1 - (P1.a)*P2'

    Ap = ap+ap*B
    Am = ap-ap*B

    assert str(Ap) == '-2*(P2.a)*P1'
    assert str(Am) == '-2*(P1.a)*P2'

    assert str(Ap*Ap) == '0'
    assert str(Am*Am) == '0'

    aB = a|B
    assert str(aB) == '-(P2.a)*P1 + (P1.a)*P2'
    return

def test_reciprocal_frame_test():

    metric = '1 # #,'+ \
             '# 1 #,'+ \
             '# # 1,'

    (e1,e2,e3) = MV.setup('e1 e2 e3',metric)


    E = e1^e2^e3
    Esq = (E*E).scalar()
    assert str(E) == 'e1^e2^e3'
    assert str(Esq) == '(e1.e2)**2 - 2*(e1.e2)*(e1.e3)*(e2.e3) + (e1.e3)**2 + (e2.e3)**2 - 1'
    Esq_inv = 1/Esq

    E1 = (e2^e3)*E
    E2 = (-1)*(e1^e3)*E
    E3 = (e1^e2)*E

    assert str(E1) == '((e2.e3)**2 - 1)*e1 + ((e1.e2) - (e1.e3)*(e2.e3))*e2 + (-(e1.e2)*(e2.e3) + (e1.e3))*e3'
    assert str(E2) == '((e1.e2) - (e1.e3)*(e2.e3))*e1 + ((e1.e3)**2 - 1)*e2 + (-(e1.e2)*(e1.e3) + (e2.e3))*e3'
    assert str(E3) == '(-(e1.e2)*(e2.e3) + (e1.e3))*e1 + (-(e1.e2)*(e1.e3) + (e2.e3))*e2 + ((e1.e2)**2 - 1)*e3'

    w = (E1|e2)
    w = w.expand()
    assert str(w) == '0'

    w = (E1|e3)
    w = w.expand()
    assert str(w) == '0'

    w = (E2|e1)
    w = w.expand()
    assert str(w) == '0'

    w = (E2|e3)
    w = w.expand()
    assert str(w) == '0'

    w = (E3|e1)
    w = w.expand()
    assert str(w) == '0'

    w = (E3|e2)
    w = w.expand()
    assert str(w) == '0'

    w = (E1|e1)
    w = (w.expand()).scalar()
    Esq = expand(Esq)
    assert str(simplify(w/Esq)) == '1'

    w = (E2|e2)
    w = (w.expand()).scalar()
    assert str(simplify(w/Esq)) == '1'

    w = (E3|e3)
    w = (w.expand()).scalar()
    assert str(simplify(w/Esq)) == '1'
    return

def test_dummy():
    return
