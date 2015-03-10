#!/usr/bin/env python

from __future__ import print_function

from sympy import symbols, log, simplify, diff, cos, sin
from sympy.galgebra import MV, ReciprocalFrame, Format
from sympy.galgebra import oprint
from sympy.galgebra import xdvi, Get_Program, Print_Function
from sympy.galgebra import Manifold


def Test_Reciprocal_Frame():
    Print_Function()
    Format()
    coords = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=coords)

    mfvar = (u, v) = symbols('u v')

    eu = ex + ey
    ev = ex - ey

    (eu_r, ev_r) = ReciprocalFrame([eu, ev])

    oprint('\\mbox{Frame}', (eu, ev), '\\mbox{Reciprocal Frame}', (eu_r, ev_r))

    print(r'%\bm{e}_{u}\cdot\bm{e}^{u} =', (eu | eu_r))
    print(r'%\bm{e}_{u}\cdot\bm{e}^{v} =', eu | ev_r)
    print(r'%\bm{e}_{v}\cdot\bm{e}^{u} =', ev | eu_r)
    print(r'%\bm{e}_{v}\cdot\bm{e}^{v} =', ev | ev_r)

    eu = ex + ey + ez
    ev = ex - ey

    (eu_r, ev_r) = ReciprocalFrame([eu, ev])

    oprint('\\mbox{Frame}', (eu, ev), '\\mbox{Reciprocal Frame}', (eu_r, ev_r))

    print(r'%\bm{e}_{u}\cdot\bm{e}^{u} =', eu | eu_r)
    print(r'%\bm{e}_{u}\cdot\bm{e}^{v} =', eu | ev_r)
    print(r'%\bm{e}_{v}\cdot\bm{e}^{u} =', ev | eu_r)
    print(r'%\bm{e}_{v}\cdot\bm{e}^{v} =', ev | ev_r)
    return


def Plot_Mobius_Strip_Manifold():
    Print_Function()
    coords = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=coords)
    mfvar = (u, v) = symbols('u v')
    X = (cos(u) + v*cos(u/2)*cos(u))*ex + (sin(u) + v*cos(u/2)*sin(u))*ey + v*sin(u/2)*ez
    MF = Manifold(X, mfvar, True, I=MV.I)
    MF.Plot2DSurface([0.0, 6.28, 48], [-0.3, 0.3, 12], surf=False, skip=[4, 4], tan=0.15)
    return


def Distorted_manifold_with_scalar_function():
    Print_Function()
    coords = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=coords)
    mfvar = (u, v) = symbols('u v')
    X = 2*u*ex + 2*v*ey + (u**3 + v**3/2)*ez
    MF = Manifold(X, mfvar, I=MV.I)

    (eu, ev) = MF.Basis()

    g = (v + 1)*log(u)
    dg = MF.Grad(g)
    print('g =', g)
    print('dg =', dg)
    print('\\eval{dg}{u=1,v=0} =', dg.subs({u: 1, v: 0}))
    G = u*eu + v*ev
    dG = MF.Grad(G)
    print('G =', G)
    print('P(G) =', MF.Proj(G))
    print('dG =', dG)
    print('P(dG) =', MF.Proj(dG))
    PS = u*v*eu ^ ev
    print('P(S) =', PS)
    print('dP(S) =', MF.Grad(PS))
    print('P(dP(S)) =', MF.Proj(MF.Grad(PS)))
    return


def Simple_manifold_with_scalar_function_derivative():
    Print_Function()
    coords = (x, y, z) = symbols('x y z')
    basis = (e1, e2, e3, grad) = MV.setup('e_1 e_2 e_3', metric='[1,1,1]', coords=coords)
    # Define surface
    mfvar = (u, v) = symbols('u v')
    X = u*e1 + v*e2 + (u**2 + v**2)*e3
    print('\\f{X}{u,v} =', X)
    MF = Manifold(X, mfvar)
    (eu, ev) = MF.Basis()
    # Define field on the surface.
    g = (v + 1)*log(u)

    print('\\f{g}{u,v} =', g)

    # Method 1: Using old Manifold routines.
    VectorDerivative = (MF.rbasis[0]/MF.E_sq)*diff(g, u) + (MF.rbasis[1]/MF.E_sq)*diff(g, v)
    print('\\eval{\\nabla g}{u=1,v=0} =', VectorDerivative.subs({u: 1, v: 0}))

    # Method 2: Using new Manifold routines.
    dg = MF.Grad(g)
    print('\\eval{\\f{Grad}{g}}{u=1,v=0} =', dg.subs({u: 1, v: 0}))
    dg = MF.grad*g
    print('\\eval{\\nabla g}{u=1,v=0} =', dg.subs({u: 1, v: 0}))
    return


def Simple_manifold_with_vector_function_derivative():
    Print_Function()
    coords = (x, y, z) = symbols('x y z')
    basis = (ex, ey, ez, grad) = \
            MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=coords)
    # Define surface
    mfvar = (u, v) = symbols('u v')
    X = u*ex + v*ey + (u**2 + v**2)*ez
    print('\\f{X}{u,v} =', X)
    MF = Manifold(X, mfvar)
    (eu, ev) = MF.Basis()

    # Define field on the surface.
    g = (v + 1)*log(u)

    print('\\mbox{Scalar Function: } g =', g)
    dg = MF.grad*g
    dg.Fmt(3, '\\mbox{Scalar Function Derivative: } \\nabla g')
    print('\\eval{\\nabla g}{(1,0)} =', dg.subs({u: 1, v: 0}))

    # Define vector field on the surface

    G = v**2*eu + u**2*ev
    print('\\mbox{Vector Function: } G =', G)
    dG = MF.grad*G
    dG.Fmt(3, '\\mbox{Vector Function Derivative: } \\nabla G')
    print('\\eval{\\nabla G}{(1,0)} =', dG.subs({u: 1, v: 0}))

    return


def dummy():
    return


def main():
    Get_Program()

    Test_Reciprocal_Frame()
    Distorted_manifold_with_scalar_function()
    Simple_manifold_with_scalar_function_derivative()
    Simple_manifold_with_vector_function_derivative()
    # Plot_Mobius_Strip_Manifold()
    xdvi()
    return

if __name__ == "__main__":
    main()
