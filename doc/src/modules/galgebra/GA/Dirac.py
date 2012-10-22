#!/usr/local/bin/python
#Dirac.py

from sympy.galgebra.GA import *
from sympy.galgebra.latex_ex import *
from sympy import *

if __name__ == '__main__':

    metric = '1  0  0  0,'+\
             '0 -1  0  0,'+\
             '0  0 -1  0,'+\
             '0  0  0 -1'

    vars = symbols('t x y z')
    gamma_t,gamma_x,gamma_y,gamma_z = MV.setup('gamma_t gamma_x gamma_y gamma_z',metric,True,vars)

    m,e = symbols('m e')
    Format('1 1 1 1')
    I = MV(ONE,'pseudo')
    nvars = len(vars)
    psi = MV('psi','spinor',fct=True)
    A = MV('A','vector',fct=True)
    sig_x = gamma_x*gamma_t
    sig_y = gamma_y*gamma_t
    sig_z = gamma_z*gamma_t
    print '$A$ is 4-vector potential'
    print A
    print r'$\bm{\psi}$ is 8-component real spinor (even multi-vector)'
    print psi
    dirac_eq = psi.grad()*I*sig_z-e*A*psi-m*psi*gamma_t
    dirac_eq.simplify()
    print 'Dirac equation in terms of real geometric algebra/calculus '+\
          r'$\lp\nabla \bm{\psi} I \sigma_{z}-eA\bm{\psi} = m\bm{\psi}\gamma_{t}\rp$'
    print 'Spin measured with respect to $z$ axis'
    Format('mv=3')
    print r'\nabla \bm{\psi} I \sigma_{z}-eA\bm{\psi}-m\bm{\psi}\gamma_{t} = ',dirac_eq,' = 0'
    xdvi(filename='Dirac.tex')
