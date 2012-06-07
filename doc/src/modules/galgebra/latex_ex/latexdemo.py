import sys, sympy
import sympy.galgebra.latex_ex as tex

if __name__ == '__main__':

    tex.Format()
    xbm,alpha_1,delta__nugamma_r = sympy.symbols('xbm alpha_1 delta__nugamma_r')

    x = alpha_1*xbm/delta__nugamma_r

    print 'x =',x

    tex.xdvi()
