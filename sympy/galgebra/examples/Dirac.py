#!/usr/bin/python
#Dirac.py

import sympy.galgebra.GA as GA
import sympy.galgebra.latex_ex as tex
import sys

GA.set_main(sys.modules[__name__])

if __name__ == '__main__':

    metric = '1  0  0  0,'+\
             '0 -1  0  0,'+\
             '0  0 -1  0,'+\
             '0  0  0 -1'

    vars = GA.make_symbols('t x y z')
    GA.MV.setup('gamma_t gamma_x gamma_y gamma_z',metric,True,vars)

    parms = GA.make_symbols('m e')
    tex.Format()
    I = GA.MV(GA.ONE,'pseudo')
    nvars = len(vars)
    psi = GA.MV('psi','spinor',fct=True)
    psi.convert_to_blades()
    A = GA.MV('A','vector',fct=True)
    sig_x = gamma_x*gamma_t
    sig_y = gamma_y*gamma_t
    sig_z = gamma_z*gamma_t
    print '$A$ is 4-vector potential'
    print A
    print r'$\bm{\psi}$ is 8-component real spinor (even multi-vector)'
    print psi
    dirac_eq = psi.grad()*I*sig_z-e*A*psi-m*psi*gamma_t
    dirac_eq.simplify()
    dirac_eq.convert_to_blades()
    print 'Dirac equation in terms of real geometric algebra/calculus '+\
          r'$\lp\nabla \bm{\psi} I \sigma_{z}-eA\bm{\psi} = m\bm{\psi}\gamma_{t}\rp$'
    print 'Spin measured with respect to $z$ axis'
    tex.MV_format(3)
    print r'\nabla \bm{\psi} I \sigma_{z}-eA\bm{\psi}-m\bm{\psi}\gamma_{t} = ',dirac_eq,' = 0'
    tex.xdvi(filename='Dirac.tex')
