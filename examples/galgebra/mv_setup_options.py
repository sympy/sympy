#!/usr/bin/env python

from __future__ import division, print_function

from sympy import symbols
from sympy.galgebra import MV
from sympy.galgebra import enhance_print, Get_Program, Print_Function

def MV_setup_options():
    Print_Function()

    (e1, e2, e3) = MV.setup('e_1 e_2 e_3', '[1,1,1]')
    v = MV('v', 'vector')
    print(v)

    (e1, e2, e3) = MV.setup('e*1|2|3', '[1,1,1]')
    v = MV('v', 'vector')
    print(v)

    (e1, e2, e3) = MV.setup('e*x|y|z', '[1,1,1]')
    v = MV('v', 'vector')
    print(v)

    coords = symbols('x y z')
    (e1, e2, e3, grad) = MV.setup('e', '[1,1,1]', coords=coords)
    v = MV('v', 'vector')
    print(v)

    return

def dummy():
    return

def main():
    Get_Program(True)
    enhance_print()
    MV_setup_options()
    return

if __name__ == "__main__":
    main()
