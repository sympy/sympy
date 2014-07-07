from sympy import sin,cos,sinh,cosh,symbols,expand,simplify
from sympy.galgebra.printer import Format,xpdf
from sympy.galgebra.mv import Com
from sympy.galgebra.ga import Ga

def main():
    Format()
    (g3d,ex,ey,ez) = Ga.build('e*x|y|z')
    A = g3d.mv('A','mv')
    print r'\bm{A} =',A
    A.Fmt(2,r'\bm{A}')
    A.Fmt(3,r'\bm{A}')

    X = (x,y,z) = symbols('x y z',real=True)
    o3d = Ga('e_x e_y e_z',g=[1,1,1],coords=X)
    (ex,ey,ez) = o3d.mv()

    f = o3d.mv('f','scalar',f=True)
    A = o3d.mv('A','vector',f=True)
    B = o3d.mv('B','bivector',f=True)

    print r'\bm{A} =',A
    print r'\bm{B} =',B

    print 'grad*f =',o3d.grad*f
    print r'grad|\bm{A} =',o3d.grad|A
    print r'grad*\bm{A} =',o3d.grad*A

    print r'-I*(grad^\bm{A}) =',-o3d.i*(o3d.grad^A)
    print r'grad*\bm{B} =',o3d.grad*B
    print r'grad^\bm{B} =',o3d.grad^B
    print r'grad|\bm{B} =',o3d.grad|B

    g4d = Ga('a b c d')

    (a,b,c,d) = g4d.mv()

    print 'g_{ij} =',g4d.g

    print '\\bm{a|(b*c)} =',a|(b*c)
    print '\\bm{a|(b^c)} =',a|(b^c)
    print '\\bm{a|(b^c^d)} =',a|(b^c^d)
    print '\\bm{a|(b^c)+c|(a^b)+b|(c^a)} =',(a|(b^c))+(c|(a^b))+(b|(c^a))
    print '\\bm{a*(b^c)-b*(a^c)+c*(a^b)} =',a*(b^c)-b*(a^c)+c*(a^b)
    print '\\bm{a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)} =',a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)
    print '\\bm{(a^b)|(c^d)} =',(a^b)|(c^d)
    print '\\bm{((a^b)|c)|d} =',((a^b)|c)|d
    print '\\bm{(a^b)\\times (c^d)} =',Com(a^b,c^d)

    g = '1 # #,'+ \
         '# 1 #,'+ \
         '# # 1'

    ng3d = Ga('e1 e2 e3',g=g)
    (e1,e2,e3) = ng3d.mv()

    E = e1^e2^e3
    Esq = (E*E).scalar()
    print 'E =',E
    print '%E^{2} =',Esq
    Esq_inv = 1/Esq

    E1 = (e2^e3)*E
    E2 = (-1)*(e1^e3)*E
    E3 = (e1^e2)*E

    print 'E1 = (e2^e3)*E =',E1
    print 'E2 =-(e1^e3)*E =',E2
    print 'E3 = (e1^e2)*E =',E3

    print 'E1|e2 =',(E1|e2).expand()
    print 'E1|e3 =',(E1|e3).expand()
    print 'E2|e1 =',(E2|e1).expand()
    print 'E2|e3 =',(E2|e3).expand()
    print 'E3|e1 =',(E3|e1).expand()
    print 'E3|e2 =',(E3|e2).expand()
    w = ((E1|e1).expand()).scalar()
    Esq = expand(Esq)
    print '%(E1\\cdot e1)/E^{2} =',simplify(w/Esq)
    w = ((E2|e2).expand()).scalar()
    print '%(E2\\cdot e2)/E^{2} =',simplify(w/Esq)
    w = ((E3|e3).expand()).scalar()
    print '%(E3\\cdot e3)/E^{2} =',simplify(w/Esq)

    X = (r,th,phi) = symbols('r theta phi')
    s3d = Ga('e_r e_theta e_phi',g=[1,r**2,r**2*sin(th)**2],coords=X,norm=True)
    (er,eth,ephi) = s3d.mv()

    f = s3d.mv('f','scalar',f=True)
    A = s3d.mv('A','vector',f=True)
    B = s3d.mv('B','bivector',f=True)

    print 'A =',A
    print 'B =',B

    print 'grad*f =',s3d.grad*f
    print 'grad|A =',s3d.grad|A
    print '-I*(grad^A) =',-s3d.i*(s3d.grad^A)
    print 'grad^B =',s3d.grad^B

    coords = symbols('t x y z')
    m4d = Ga('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)
    (g0,g1,g2,g3) = m4d.mv()
    I = m4d.i

    B = m4d.mv('B','vector',f=True)
    E = m4d.mv('E','vector',f=True)
    B.set_coef(1,0,0)
    E.set_coef(1,0,0)
    B *= g0
    E *= g0
    J = m4d.mv('J','vector',f=True)
    F = E+I*B

    print 'B = \\bm{B\\gamma_{t}} =',B
    print 'E = \\bm{E\\gamma_{t}} =',E
    print 'F = E+IB =',F
    print 'J =',J
    gradF = m4d.grad*F
    gradF.Fmt(3,'grad*F')

    print 'grad*F = J'
    (gradF.get_grade(1)-J).Fmt(3,'%\\grade{\\nabla F}_{1} -J = 0')
    (gradF.get_grade(3)).Fmt(3,'%\\grade{\\nabla F}_{3} = 0')

    (alpha,beta,gamma) = symbols('alpha beta gamma')

    (x,t,xp,tp) = symbols("x t x' t'")
    m2d = Ga('gamma*t|x',g=[1,-1])
    (g0,g1) = m2d.mv()

    R = cosh(alpha/2)+sinh(alpha/2)*(g0^g1)
    X = t*g0+x*g1
    Xp = tp*g0+xp*g1
    print 'R =',R

    print r"#%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} = t'\bm{\gamma'_{t}}+x'\bm{\gamma'_{x}} = R\lp t'\bm{\gamma_{t}}+x'\bm{\gamma_{x}}\rp R^{\dagger}"

    Xpp = R*Xp*R.rev()
    Xpp = Xpp.collect()
    Xpp = Xpp.trigsimp()
    print r"%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} =",Xpp
    Xpp = Xpp.subs({sinh(alpha):gamma*beta,cosh(alpha):gamma})

    print r'%\f{\sinh}{\alpha} = \gamma\beta'
    print r'%\f{\cosh}{\alpha} = \gamma'

    print r"%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} =",Xpp.collect()

    coords = symbols('t x y z')
    m4d = Ga('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)
    (g0,g1,g2,g3) = m4d.mv()
    I = m4d.i
    (m,e) = symbols('m e')

    psi = m4d.mv('psi','spinor',f=True)
    A = m4d.mv('A','vector',f=True)
    sig_z = g3*g0
    print '\\bm{A} =',A
    print '\\bm{\\psi} =',psi

    dirac_eq = (m4d.grad*psi)*I*sig_z-e*A*psi-m*psi*g0
    dirac_eq.simplify()

    dirac_eq.Fmt(3,r'\nabla \bm{\psi} I \sigma_{z}-e\bm{A}\bm{\psi}-m\bm{\psi}\gamma_{t} = 0')

    xpdf()
    return

if __name__ == "__main__":
    main()
