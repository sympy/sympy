from GA import *
from GAPrint import Get_Program,Print_Function

enhance_print()

def MV_setup_options():
    Print_Function()

    (e1,e2,e3) = MV.setup('e_1 e_2 e_3','[1,1,1]')
    v = MV('v', 'vector')
    print v

    (e1,e2,e3) = MV.setup('e*1|2|3','[1,1,1]')
    v = MV('v', 'vector')
    print v

    (e1,e2,e3) = MV.setup('e*x|y|z','[1,1,1]')
    v = MV('v', 'vector')
    print v

    coords = symbols('x y z')
    (e1,e2,e3,grad) = MV.setup('e','[1,1,1]',coords=coords)
    v = MV('v', 'vector')
    print v

    return

def dummy():
    return

Get_Program()
MV_setup_options()
