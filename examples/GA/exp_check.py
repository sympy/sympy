import sys

from sympy import symbols,sin,cos
from sympy.GA.GA import MV
from sympy.GA.GAPrint import enhance_print
enhance_print()

def exp_operation():
    (ex,ey,ez) = MV.setup('e*x|y|z',metric='[1,1,1]')

    u = MV('u','vector')
    v = MV('v','vector')
    w = MV('w','vector')
    print u
    print v

    uv = u^v
    print uv
    print uv.is_blade()

    exp_uv = uv.exp()
    exp_uv.Fmt(2,'exp(uv)')

    return

exp_operation()
