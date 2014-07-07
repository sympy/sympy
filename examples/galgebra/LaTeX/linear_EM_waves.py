import sys
from sympy import symbols,sin,cos,exp,I,Matrix,solve,simplify
from sympy.galgebra.printer import Format,xpdf,Get_Program,Print_Function
from sympy.galgebra.ga import Ga
from sympy.galgebra.metric import linear_expand

def EM_Waves_in_Geom_Calculus_Complex():
    #Print_Function()
    X = (t,x,y,z) = symbols('t x y z',real=True)
    g = '1 # # 0,# 1 # 0,# # 1 0,0 0 0 -1'
    coords = (xE,xB,xk,t) = symbols('x_E x_B x_k t',real=True)
    (EBkst,eE,eB,ek,et) = Ga.build('e_E e_B e_k e_t',g=g,coords=coords)

    i = EBkst.i

    E,B,k,w = symbols('E B k omega',real=True)

    F = E*eE*et+i*B*eB*et
    K = k*ek+w*et
    X = xE*eE+xB*eB+xk*ek+t*et
    KX = (K|X).scalar()
    F = F*exp(I*KX)

    g = EBkst.g

    print 'g =', g
    print 'X =', X
    print 'K =', K
    print 'K|X =', KX
    print 'F =', F

    gradF = EBkst.grad*F

    gradF = gradF.simplify()

    (gradF).Fmt(3,'grad*F = 0')

    gradF = gradF.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    KX = KX.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    print r'%\mbox{Substituting }e_{E}\cdot e_{B} = e_{E}\cdot e_{k} = e_{B}\cdot e_{k} = 0'

    (gradF / (I*exp(I*KX))).Fmt(3,r'%\lp\bm{\nabla}F\rp/\lp ie^{iK\cdot X}\rp = 0')

    return

def EM_Waves_in_Geom_Calculus_Real():
    #Print_Function()
    X = (t,x,y,z) = symbols('t x y z',real=True)
    g = '1 # # 0,# 1 # 0,# # 1 0,0 0 0 -1'
    coords = (xE,xB,xk,t) = symbols('x_E x_B x_k t',real=True)
    (EBkst,eE,eB,ek,et) = Ga.build('e_E e_B e_k e_t',g=g,coords=coords)

    i = EBkst.i

    E,B,k,w = symbols('E B k omega',real=True)

    F = E*eE*et+i*B*eB*et
    K = k*ek+w*et
    X = xE*eE+xB*eB+xk*ek+t*et
    KX = (K|X).scalar()
    F = F*sin(KX)

    g = EBkst.g

    print 'g =', g
    print 'X =', X
    print 'K =', K
    print 'K|X =', KX
    F.Fmt(3,'F')

    gradF = EBkst.grad*F

    gradF = gradF.simplify()

    (gradF).Fmt(3,'grad*F = 0')

    gradF = gradF.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    KX = KX.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    print r'%\mbox{Substituting }e_{E}\cdot e_{B} = e_{E}\cdot e_{k} = e_{B}\cdot e_{k} = 0'

    (gradF / (cos(KX))).Fmt(3,r'%\lp\bm{\nabla}F\rp/\lp \f{\cos}{K\cdot X}\rp = 0')

    return

def dummy():
    return

def main():
    #Get_Program()
    Format()

    #EM_Waves_in_Geom_Calculus_Complex()
    EM_Waves_in_Geom_Calculus_Real()
    xpdf()
    return

if __name__ == "__main__":
    main()
