#!/usr/bin/python

import sys
import sympy.galgebra.GA as GA
import sympy.galgebra.latex_ex as tex

GA.set_main(sys.modules[__name__])

if __name__ == '__main__':

    metric = '1  0  0  0,'+\
             '0 -1  0  0,'+\
             '0  0 -1  0,'+\
             '0  0  0 -1'

    vars = GA.make_symbols('t x y z')
    GA.MV.setup('gamma_t gamma_x gamma_y gamma_z',metric,True,vars)
    tex.Format()
    I = GA.MV(1,'pseudo')
    I.convert_to_blades()
    print '$I$ Pseudo-Scalar'
    print 'I =',I
    B = GA.MV('B','vector',fct=True)
    E = GA.MV('E','vector',fct=True)
    B.set_coef(1,0,0)
    E.set_coef(1,0,0)
    B *= gamma_t
    E *= gamma_t
    B.convert_to_blades()
    E.convert_to_blades()
    J = GA.MV('J','vector',fct=True)
    print '$B$ Magnetic Field Bi-Vector'
    print 'B = Bvec gamma_0 =',B
    print '$E$ Electric Field Bi-Vector'
    print 'E = Evec gamma_0 =',E
    F = E+I*B
    print '$E+IB$ Electo-Magnetic Field Bi-Vector'
    print 'F = E+IB =',F
    print '$J$ Four Current'
    print 'J =',J
    gradF = F.grad()
    gradF.convert_to_blades()
    print 'Geometric Derivative of Electo-Magnetic Field Bi-Vector'
    tex.MV_format(3)
    print '\\nabla F =',gradF
    print 'All Maxwell Equations are'
    print '\\nabla F = J'
    print 'Div $E$ and Curl $H$ Equations'
    print '<\\nabla F>_1 -J =',gradF.project(1)-J,' = 0'
    print 'Curl $E$ and Div $B$ equations'
    print '<\\nabla F>_3 =',gradF.project(3),' = 0'
    tex.xdvi(filename='Maxwell.tex')
