import sys

from sympy import symbols,sin,cos
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Eprint

def main():
    Eprint()
    (o3d,ex,ey,ez) = Ga.build('e*x|y|z',g=[1,1,1])

    u = o3d.mv('u','vector')
    v = o3d.mv('v','vector')
    w = o3d.mv('w','vector')
    print u
    print v

    uv = u^v
    print uv
    print uv.is_blade()

    exp_uv = uv.exp()
    print 'exp(uv) =', exp_uv

    return

if __name__ == "__main__":
    main()
