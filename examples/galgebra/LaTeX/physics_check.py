
from sympy import symbols,sin,cos
from sympy.galgebra.printer import Format,xpdf,Get_Program,Print_Function
from sympy.galgebra.ga import Ga

def Maxwells_Equations_in_Geom_Calculus():
    Print_Function()
    X = symbols('t x y z',real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=X)

    I = st4d.i

    B = st4d.mv('B','vector',f=True)
    E = st4d.mv('E','vector',f=True)
    B.set_coef(1,0,0)
    E.set_coef(1,0,0)
    B *= g0
    E *= g0
    J = st4d.mv('J','vector',f=True)
    F = E+I*B

    print r'\text{Pseudo Scalar\;\;}I =',I
    print '\\text{Magnetic Field Bi-Vector\\;\\;} B = \\bm{B\\gamma_{t}} =',B
    print '\\text{Electric Field Bi-Vector\\;\\;} E = \\bm{E\\gamma_{t}} =',E
    print '\\text{Electromagnetic Field Bi-Vector\\;\\;} F = E+IB =',F
    print '%\\text{Four Current Density\\;\\;} J =',J
    gradF = st4d.grad*F
    print '#Geom Derivative of Electomagnetic Field Bi-Vector'
    gradF.Fmt(3,'grad*F')

    print '#Maxwell Equations'
    print 'grad*F = J'
    print '#Div $E$ and Curl $H$ Equations'
    (gradF.get_grade(1)-J).Fmt(3,'%\\grade{\\nabla F}_{1} -J = 0')
    print '#Curl $E$ and Div $B$ equations'
    (gradF.get_grade(3)).Fmt(3,'%\\grade{\\nabla F}_{3} = 0')
    return

def Dirac_Equation_in_Geog_Calculus():
    Print_Function()
    coords = symbols('t x y z',real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)
    I = st4d.i

    (m,e) = symbols('m e')

    psi = st4d.mv('psi','spinor',f=True)
    A = st4d.mv('A','vector',f=True)
    sig_z = g3*g0

    print '\\text{4-Vector Potential\\;\\;}\\bm{A} =',A
    print '\\text{8-component real spinor\\;\\;}\\bm{\\psi} =',psi

    dirac_eq = (st4d.grad*psi)*I*sig_z-e*A*psi-m*psi*g0
    dirac_eq = dirac_eq.simplify()

    dirac_eq.Fmt(3,r'%\text{Dirac Equation\;\;}\nabla \bm{\psi} I \sigma_{z}-e\bm{A}\bm{\psi}-m\bm{\psi}\gamma_{t} = 0')

    return

def Lorentz_Tranformation_in_Geog_Algebra():
    Print_Function()
    (alpha,beta,gamma) = symbols('alpha beta gamma')
    (x,t,xp,tp) = symbols("x t x' t'",real=True)
    (st2d,g0,g1) = Ga.build('gamma*t|x',g=[1,-1])

    from sympy import sinh,cosh

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
    return

def dummy():
    return

def main():
    Get_Program()
    Format()

    Maxwells_Equations_in_Geom_Calculus()
    Dirac_Equation_in_Geog_Calculus()
    Lorentz_Tranformation_in_Geog_Algebra()

    xpdf()
    return

if __name__ == "__main__":
    main()
