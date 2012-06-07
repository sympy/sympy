
from sympy import *
from sympy.galgebra.GA import *
from sympy.galgebra.latex_ex import *

if __name__ == '__main__':

    metric = '1  0  0  0,'+\
             '0 -1  0  0,'+\
             '0  0 -1  0,'+\
             '0  0  0 -1'

    vars = symbols('t x y z')
    gamma_t,gamma_x,gamma_y,gamma_z = MV.setup('gamma_t gamma_x gamma_y gamma_z',metric,True,vars)
    LatexPrinter.format(1,1,1,1)
    I = MV(1,'pseudo')
    print '$I$ Pseudo-Scalar'
    print 'I =',I
    B = MV('B','vector',fct=True)
    E = MV('E','vector',fct=True)
    B.set_coef(1,0,0)
    E.set_coef(1,0,0)
    B *= gamma_t
    E *= gamma_t
    J = MV('J','vector',fct=True)
    F = E+I*B
    print ' '
    print '$B$ Magnetic Field Bi-Vector'
    print 'B = Bvec gamma_0 =',B
    print '$F$ Electric Field Bi-Vector'
    print 'E = Evec gamma_0 =',E
    print '$E+IB$ Electo-Magnetic Field Bi-Vector'
    print 'F = E+IB =',F
    print '$J$ Four Current'
    print 'J =',J
    gradF = F.grad()
    print 'Geometric Derivative of Electo-Magnetic Field Bi-Vector'
    MV_format(3)
    print '\\nabla F =',gradF
    print 'All Maxwell Equations are'
    print '\\nabla F = J'
    print 'Div $E$ and Curl $H$ Equations'
    print '<\\nabla F>_1 -J =',gradF.project(1)-J,' = 0'
    print 'Curl $E$ and Div $B$ equations'
    print '<\\nabla F>_3 =',gradF.project(3),' = 0'
    xdvi(filename='Maxwell.tex')
