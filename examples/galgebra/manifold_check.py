#!/usr/bin/env python

from __future__ import division, print_function

from sympy import symbols, log, simplify, diff, cos, sin
from sympy.galgebra.ga import MV, ReciprocalFrame
from sympy.galgebra.debug import oprint
from sympy.galgebra.printing import GA_Printer, enhance_print, Get_Program, Print_Function
from sympy.galgebra.manifold import Manifold


def Test_Reciprocal_Frame():
    Print_Function()
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
    return


def Plot_Mobius_Strip_Manifold():
    Print_Function()
    coords = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('ex ey ez', metric='[1,1,1]', coords=coords)
    mfvar = (u, v) = symbols('u v')
    X = (cos(u) + v*cos(u/2)*cos(u))*ex + (sin(u) + v*cos(u/2)*sin(u))*ey + v*sin(u/2)*ez
    MF = Manifold(X, mfvar, True, I=MV.I)
    MF.Plot2DSurface([0.0, 6.28, 48], [-0.3, 0.3, 12], surf=False, skip=[4, 4], tan=0.15)
    return


def Distorted_manifold_with_scalar_function():
    Print_Function()
    coords = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('ex ey ez', metric='[1,1,1]', coords=coords)
    mfvar = (u, v) = symbols('u v')
    X = 2*u*ex + 2*v*ey + (u**3 + v**3/2)*ez
    MF = Manifold(X, mfvar, I=MV.I)

    (eu, ev) = MF.Basis()

    g = (v + 1)*log(u)
    dg = MF.Grad(g)
    print('g =', g)
    print('dg =', dg)
    print('dg(1,0) =', dg.subs({u: 1, v: 0}))
    G = u*eu + v*ev
    dG = MF.Grad(G)
    print('G =', G)
    print('P(G) =', MF.Proj(G))
    print('zcoef =', simplify(2*(u**2 + v**2)*(-4*u**2 - 4*v**2 - 1)))
    print('dG =', dG)
    print('P(dG) =', MF.Proj(dG))
    PS = u*v*eu ^ ev
    print('PS =', PS)
    print('dPS =', MF.Grad(PS))
    print('P(dPS) =', MF.Proj(MF.Grad(PS)))
    return


def Simple_manifold_with_scalar_function_derivative():
    Print_Function()
    coords = (x, y, z) = symbols('x y z')
    basis = (e1, e2, e3, grad) = MV.setup('e_1 e_2 e_3', metric='[1,1,1]', coords=coords)
    # Define surface
    mfvar = (u, v) = symbols('u v')
    X = u*e1 + v*e2 + (u**2 + v**2)*e3
    print(X)
    MF = Manifold(X, mfvar)

    # Define field on the surface.
    g = (v + 1)*log(u)

    # Method 1: Using old Manifold routines.
    VectorDerivative = (MF.rbasis[0]/MF.E_sq)*diff(g, u) + (MF.rbasis[1]/MF.E_sq)*diff(g, v)
    print('Vector derivative =', VectorDerivative.subs({u: 1, v: 0}))

    # Method 2: Using new Manifold routines.
    dg = MF.Grad(g)
    print('Vector derivative =', dg.subs({u: 1, v: 0}))
    return


def dummy():
    return


def main():
    Get_Program(True)
    with GA_Printer():
        enhance_print()
        Test_Reciprocal_Frame()
        Distorted_manifold_with_scalar_function()
        Simple_manifold_with_scalar_function_derivative()
        # Plot_Mobius_Strip_Manifold()
    return

if __name__ == "__main__":
    main()
