#!/usrlocal/bin/python
#EandM.py

from sympy.galgebra.GA import *
from sympy.galgebra.latex_ex import *
from sympy import *

import sympy,numpy,sys

if __name__ == '__main__':
    metric = '1 0 0,'+\
             '0 1 0,'+\
             '0 0 1'

    gamma_x,gamma_y,gamma_z = MV.setup('gamma_x gamma_y gamma_z',metric,True)
    Format('1 1 1 1')

    coords = r,theta,phi = symbols('r theta phi')
    x = r*(sympy.cos(theta)*gamma_z+sympy.sin(theta)*\
        (sympy.cos(phi)*gamma_x+sympy.sin(phi)*gamma_y))
    x.set_name('x')

    MV.rebase(x,coords,'e',False)

    #psi = MV.scalar_fct('psi')
    psi = MV('psi','scalar',fct=True)
    #psi.name = 'psi'
    dpsi = psi.grad()
    print 'Gradient of Scalar Function $\\psi$'
    print '\\nabla\\psi =',dpsi

    #A = MV.vector_fct('A')
    A = MV('A','vector',fct=True)
    #A.name = 'A'
    print 'Div and Curl of Vector Function $A$'
    print A

    gradA = A.grad()
    I = MV(ONE,'pseudo')
    divA = A.grad_int()
    curlA = -I*A.grad_ext()
    print '\\nabla \\cdot A =',divA
    Format('mv=3')
    print '-I\\lp\\nabla \\W A\\rp =',curlA

    xdvi(filename='coords.tex')
