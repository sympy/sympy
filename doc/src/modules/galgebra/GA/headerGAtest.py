import os,sys,sympy
from sympy.galgebra.GA import set_main, make_symbols, types, MV, ZERO, ONE, HALF
from sympy import collect
set_main(sys.modules[__name__])

def F(x):
        """
        Conformal Mapping Function
        """
        Fx = HALF*((x*x)*n+2*x-nbar)
        return(Fx)

def make_vector(a,n = 3):
        if type(a) == types.StringType:
                sym_str = ''
                for i in range(n):
                        sym_str += a+str(i)+' '
                sym_lst = make_symbols(sym_str)
                sym_lst.append(ZERO)
                sym_lst.append(ZERO)
                a = MV(sym_lst,'vector')
        return(F(a))

if __name__ == '__main__':
