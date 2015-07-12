#!/usr/bin/env python

from __future__ import division, print_function

from sympy import symbols, sin, cos
from sympy.galgebra import MV
from sympy.galgebra import enhance_print

def main():
    enhance_print()
    (ex, ey, ez) = MV.setup('e*x|y|z', metric='[1,1,1]')

    u = MV('u', 'vector')
    v = MV('v', 'vector')
    w = MV('w', 'vector')
    print(u)
    print(v)

    uv = u ^ v
    print(uv)
    print(uv.is_blade())

    exp_uv = uv.exp()
    exp_uv.Fmt(2, 'exp(uv)')

    return

if __name__ == "__main__":
    main()
