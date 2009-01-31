import sys
import sympy.galgebra.GA as GA
import sympy.galgebra.latex_ex as tex

GA.set_main(sys.modules[__name__])

if __name__ == '__main__':

    tex.Format()
    GA.make_symbols('xbm alpha_1 delta__nugamma_r')
    
    x = alpha_1*xbm/delta__nugamma_r
    
    print 'x =',x
    
    tex.xdvi()
