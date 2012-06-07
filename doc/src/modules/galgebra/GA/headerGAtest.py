import os,sys,sympy
from sympy.galgebra.GA import MV, ZERO, ONE, HALF
from sympy import collect, symbols

def F(x, n, nbar):
        """
        Conformal Mapping Function
        """
        Fx = HALF*((x*x)*n+2*x-nbar)
        return(Fx)

if __name__ == '__main__':
