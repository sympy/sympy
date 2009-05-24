#!/usr/bin/python
#EandM.py

import sys
import sympy.galgebra.GA as GA
import sympy.galgebra.latex_ex as tex
import sympy,numpy,time

GA.set_main(sys.modules[__name__])

if __name__ == '__main__':
    metric = '1 0 0,'+\
             '0 1 0,'+\
             '0 0 1'

    GA.MV.setup('gamma_x gamma_y gamma_z',metric,True)
    tex.Format()

    coords = GA.make_symbols('r theta phi')
    x = r*(sympy.cos(theta)*gamma_z+sympy.sin(theta)*\
        (sympy.cos(phi)*gamma_x+sympy.sin(phi)*gamma_y))
    x.set_name('x')

    GA.MV.rebase(x,coords,'e',True)

    psi = GA.MV('psi','scalar',fct=True)

    dpsi = psi.grad()
    print 'Gradient of Scalar Function $\\psi$'
    print '\\nabla\\psi =',dpsi

    A = GA.MV('A','vector',fct=True)

    print 'Div and Curl of Vector Function $A$'
    print A

    gradA = A.grad()
    I = GA.MV(GA.ONE,'pseudo')
    divA = A.grad_int()
    curlA = -I*A.grad_ext()
    print '\\nabla \\cdot A =',divA
    tex.MV_format(3)
    print '-I\\lp\\nabla \\W A\\rp =',curlA

    tex.xdvi(filename='coords.tex')
