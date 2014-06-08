#!/usr/bin/env python

from __future__ import print_function

from sympy import symbols
from sympy.galgebra import MV, Format
from sympy.galgebra import xdvi

def main():
    Format()

    coords = (x, y, z) = symbols('x y z')

    (ex, ey, ez, grad) = MV.setup('e*x|y|z', '[1,1,1]', coords=coords)

    s = MV('s', 'scalar')
    v = MV('v', 'vector')
    b = MV('b', 'bivector')

    print(r'#3D Orthogonal Metric\newline')

    print('#Multvectors:')
    print('s =', s)
    print('v =', v)
    print('b =', b)

    print('#Products:')

    X = ((s, 's'), (v, 'v'), (b, 'b'))

    for xi in X:
        print('')
        for yi in X:
            print(xi[1] + '*' + yi[1] + ' =', xi[0]*yi[0])
            print(xi[1] + '^' + yi[1] + ' =', xi[0] ^ yi[0])
            print(xi[1] + '|' + yi[1] + ' =', xi[0] | yi[0])
            print(xi[1] + '<' + yi[1] + ' =', xi[0] < yi[0])
            print(xi[1] + '>' + yi[1] + ' =', xi[0] > yi[0])

    fs = MV('s', 'scalar', fct=True)
    fv = MV('v', 'vector', fct=True)
    fb = MV('b', 'bivector', fct=True)

    print('#Multivector Functions:')

    print('s(X) =', fs)
    print('v(X) =', fv)
    print('b(X) =', fb)

    print('#Products:')

    fX = ((grad, 'grad'), (fs, 's'), (fv, 'v'), (fb, 'b'))

    for xi in fX:
        print('')
        for yi in fX:
            if xi[1] == 'grad' and yi[1] == 'grad':
                pass
            else:
                print(xi[1] + '*' + yi[1] + ' =', xi[0]*yi[0])
                print(xi[1] + '^' + yi[1] + ' =', xi[0] ^ yi[0])
                print(xi[1] + '|' + yi[1] + ' =', xi[0] | yi[0])
                print(xi[1] + '<' + yi[1] + ' =', xi[0] < yi[0])
                print(xi[1] + '>' + yi[1] + ' =', xi[0] > yi[0])

    (ex, ey, grad) = MV.setup('e', coords=(x, y))

    print(r'#General 2D Metric\newline')
    print('#Multivector Functions:')

    s = MV('s', 'scalar', fct=True)
    v = MV('v', 'vector', fct=True)
    b = MV('v', 'bivector', fct=True)

    print('s(X) =', s)
    print('v(X) =', v)
    print('b(X) =', b)

    X = ((grad, 'grad'), (s, 's'), (v, 'v'))

    print('#Products:')

    for xi in X:
        print('')
        for yi in X:
            if xi[1] == 'grad' and yi[1] == 'grad':
                pass
            else:
                print(xi[1] + '*' + yi[1] + ' =', xi[0]*yi[0])
                print(xi[1] + '^' + yi[1] + ' =', xi[0] ^ yi[0])
                print(xi[1] + '|' + yi[1] + ' =', xi[0] | yi[0])
                print(xi[1] + '<' + yi[1] + ' =', xi[0] < yi[0])
                print(xi[1] + '>' + yi[1] + ' =', xi[0] > yi[0])

    xdvi(paper='letter')
    return

if __name__ == "__main__":
    main()
