#!/usr/bin/python

import sys
import os,sympy,time
from sympy.galgebra.GA import set_main, make_symbols, types, MV, ZERO, ONE, HALF
from sympy import collect
set_main(sys.modules[__name__])

def F(x):
    """
    Conformal Mapping Function
    """
    Fx = HALF*((x*x)*n+2*x-nbar)
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

if __name__ == '__main__':

    ti = time.time()

    MV.setup('a b c d e')
    MV.set_str_format(1)

    print 'e|(a^b) =',e|(a^b)
    print 'e|(a^b^c) =',e|(a^b^c)
    print 'a*(b^c)-b*(a^c)+c*(a^b) =',a*(b^c)-b*(a^c)+c*(a^b)
    print 'e|(a^b^c^d) =',e|(a^b^c^d)
    print -d*(a^b^c)+c*(a^b^d)-b*(a^c^d)+a*(b^c^d)

    print (a^b)|(c^d)

    # FIXME: currently broken
    """
    print 'Example: non-euclidian distance calculation'

    metric = '0 # #,# 0 #,# # 1'
    MV.setup('X Y e',metric)
    MV.set_str_format(1)
    L = X^Y^e
    B = L*e
    Bsq = (B*B)()
    print 'L = X^Y^e is a non-euclidian line'
    print 'B = L*e =',B
    BeBr =B*e*B.rev()
    print 'B*e*B.rev() =',BeBr
    print 'B^2 =',Bsq
    print 'L^2 =',(L*L)()
    make_symbols('s c Binv M S C alpha')
    Bhat = Binv*B # Normalize translation generator
    R = c+s*Bhat # Rotor R = exp(alpha*Bhat/2)
    print 's = sinh(alpha/2) and c = cosh(alpha/2)'
    print 'R = exp(alpha*B/(2*|B|)) =',R
    Z = R*X*R.rev()
    Z.expand()
    Z.collect([Binv,s,c,XdotY])
    print 'R*X*R.rev() =',Z
    W = Z|Y
    W.expand()
    W.collect([s*Binv])
    print '(R*X*rev(R)).Y =',W
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
    print '(R*X*R.rev()).Y =',W
    nl = '\n'

    Wd = collect(W,[C,S],exact=True,evaluate=False)
    print 'Wd =',Wd
    Wd_1 = Wd[ONE]
    Wd_C = Wd[C]
    Wd_S = Wd[S]
    print '|B| =',Bmag
    Wd_1 = Wd_1.subs(Bmag,1/Binv)
    Wd_C = Wd_C.subs(Bmag,1/Binv)
    Wd_S = Wd_S.subs(Bmag,1/Binv)
    print 'Wd[ONE] =',Wd_1
    print 'Wd[C] =',Wd_C
    print 'Wd[S] =',Wd_S

    lhs = Wd_1+Wd_C*C
    rhs = -Wd_S*S
    lhs = lhs**2
    rhs = rhs**2
    W = (lhs-rhs).expand()
    W = (W.subs(1/Binv**2,Bmag**2)).expand()
    print 'W =',W
    W = (W.subs(S**2,C**2-1)).expand()
    print 'W =',W
    W = collect(W,[C,C**2],evaluate=False)
    print 'W =',W

    a = W[C**2]
    b = W[C]
    c = W[ONE]

    print 'a =',a
    print 'b =',b
    print 'c =',c

    D = (b**2-4*a*c).expand()
    print 'Setting to 0 and solving for C gives:'
    print 'Discriminant D = b^2-4*a*c =',D
    C = (-b/(2*a)).expand()
    print 'C = cosh(alpha) = -b/(2*a) =',C
    """
    print '\nExample: Conformal representations of circles, lines, spheres, and planes'

    metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

    MV.setup('e0 e1 e2 n nbar',metric,debug=0)
    MV.set_str_format(1)
    e = n+nbar
    #conformal representation of points

    A = make_vector(e0)    # point a = (1,0,0)  A = F(a)
    B = make_vector(e1)    # point b = (0,1,0)  B = F(b)
    C = make_vector(-1*e0) # point c = (-1,0,0) C = F(c)
    D = make_vector(e2)    # point d = (0,0,1)  D = F(d)
    X = make_vector('x',3)

    print 'a = e0, b = e1, c = -e0, and d = e2'
    print 'A = F(a) = 1/2*(a*a*n+2*a-nbar), etc.'
    print 'Circle through a, b, and c'
    print 'Circle: A^B^C^X = 0 =',(A^B^C^X)
    print 'Line through a and b'
    print 'Line  : A^B^n^X = 0 =',(A^B^n^X)
    print 'Sphere through a, b, c, and d'
    print 'Sphere: A^B^C^D^X = 0 =',(A^B^C^D^X)
    print 'Plane through a, b, and d'
    print 'Plane : A^B^n^D^X = 0 =',(A^B^n^D^X)

    L = (A^B^e)^X

    print 'Hyperbolic Circle: (A^B^e)^X = 0 =',L

    #MV.LaTeX()

    metric = '# # # 0 0,'+ \
                     '# # # 0 0,'+ \
                     '# # # 0 0,'+ \
                     '0 0 0 0 2,'+ \
                     '0 0 0 2 0'

    MV.setup('p1 p2 p3 n nbar',metric,debug=0)
    MV.set_str_format(1)

    P1 = F(p1)
    P2 = F(p2)
    P3 = F(p3)

    print '\nExtracting direction of line from L = P1^P2^n'

    L = P1^P2^n
    delta = (L|n)|nbar
    print '(L.n).nbar=',delta

    print '\nExtracting plane of circle from C = P1^P2^P3'

    C = P1^P2^P3
    delta = ((C^n)|n)|nbar
    print '((C^n).n).nbar=',delta
    print '(p2-p1)^(p3-p1)=',(p2-p1)^(p3-p1)

    metric = '1 # #,'+ \
                     '# 1 #,'+ \
                     '# # 1,'

    MV.setup('e1 e2 e3',metric)

    print 'Example: Reciprocal Frames e1, e2, and e3 unit vectors.\n\n'

    E = e1^e2^e3
    Esq = (E*E)()
    print 'E =',E
    print 'E^2 =',Esq
    Esq_inv = 1/Esq

    E1 = (e2^e3)*E
    E2 = (-1)*(e1^e3)*E
    E3 = (e1^e2)*E

    print 'E1 = (e2^e3)*E =',E1
    print 'E2 =-(e1^e3)*E =',E2
    print 'E3 = (e1^e2)*E =',E3

    w = (E1|e2)
    w.collect(MV.g)
    w = w().expand()
    print 'E1|e2 =',w

    w = (E1|e3)
    w.collect(MV.g)
    w = w().expand()
    print 'E1|e3 =',w

    w = (E2|e1)
    w.collect(MV.g)
    w = w().expand()
    print 'E2|e1 =',w

    w = (E2|e3)
    w.collect(MV.g)
    w = w().expand()
    print 'E2|e3 =',w

    w = (E3|e1)
    w.collect(MV.g)
    w = w().expand()
    print 'E3|e1 =',w

    w = (E3|e2)
    w.collect(MV.g)
    w = w().expand()
    print 'E3|e2 =',w

    w = (E1|e1)
    w = w().expand()
    Esq = Esq.expand()
    print '(E1|e1)/E^2 =',w/Esq

    w = (E2|e2)
    w = w().expand()
    print '(E2|e2)/E^2 =',w/Esq

    w = (E3|e3)
    w = w().expand()
    print '(E3|e3)/E^2 =',w/Esq

    print '\nExtracting vectors from conformal 2 blade B = P1^P2'

    metric = ' 0 -1 #,'+ \
                     '-1  0 #,'+ \
                     ' #  # #,'

    MV.setup('P1 P2 a',metric)

    B = P1^P2
    Bsq = B*B
    print 'B^2 =',Bsq
    ap = a-(a^B)*B
    print "a' = a-(a^B)*B =",ap

    Ap = ap+ap*B
    Am = ap-ap*B

    print "A+ = a'+a'*B =",Ap
    print "A- = a'-a'*B =",Am

    print '(A+)^2 =',Ap*Ap
    print '(A-)^2 =',Am*Am

    aB = a|B
    print 'a.B =',aB

    tf = time.time()

    print 1000.0*(tf-ti)
