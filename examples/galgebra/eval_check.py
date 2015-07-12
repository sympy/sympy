#!/usr/bin/env python

from __future__ import division, print_function

from sympy import symbols
from sympy.galgebra import MV, ReciprocalFrame
from sympy.galgebra import enhance_print
from sympy.galgebra import oprint
from sympy.galgebra import define_precedence, GAeval

def main():
    enhance_print()

    coords = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('ex ey ez', metric='[1,1,1]', coords=coords)

    mfvar = (u, v) = symbols('u v')

    eu = ex + ey
    ev = ex - ey

    (eu_r, ev_r) = ReciprocalFrame([eu, ev])

    oprint('Frame', (eu, ev), 'Reciprocal Frame', (eu_r, ev_r))

    print('eu.eu_r =', eu | eu_r)
    print('eu.ev_r =', eu | ev_r)
    print('ev.eu_r =', ev | eu_r)
    print('ev.ev_r =', ev | ev_r)

    eu = ex + ey + ez
    ev = ex - ey

    (eu_r, ev_r) = ReciprocalFrame([eu, ev])

    oprint('Frame', (eu, ev), 'Reciprocal Frame', (eu_r, ev_r))

    print('eu.eu_r =', eu | eu_r)
    print('eu.ev_r =', eu | ev_r)
    print('ev.eu_r =', ev | eu_r)
    print('ev.ev_r =', ev | ev_r)

    print('eu =', eu)
    print('ev =', ev)

    define_precedence(locals())

    print(GAeval('eu^ev|ex', True))
    print(GAeval('eu^ev|ex*eu', True))
    return

if __name__ == "__main__":
    main()
