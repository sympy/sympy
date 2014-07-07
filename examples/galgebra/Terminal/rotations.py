import sys

from sympy import symbols,sin,cos,sqrt
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Eprint

def main():
    Eprint()
    (o3d,ex,ey,ez) = Ga.build('e*x|y|z',g=[1,1,1])

    (r,th,phi,alpha,beta,gamma) = symbols('r theta phi alpha beta gamma',real=True)
    (x_a,y_a,z_a,x_b,y_b,z_b,ab_mag,th_ab) = symbols('x_a y_a z_a x_b y_b z_b ab_mag theta_ab',real=True)

    I = ex^ey^ez

    a = o3d.mv('a','vector')
    b = o3d.mv('b','vector')
    c = o3d.mv('c','vector')
    ab = a-b

    print 'a =', a
    print 'b =', b
    print 'c =', c
    print 'ab =', ab

    ab_norm = ab/ab_mag

    print 'ab/|ab| =', ab_norm

    R_ab = cos(th_ab/2) +I*ab_norm*cos(th_ab/2)
    R_ab_rev = R_ab.rev()
    print 'R_ab =', R_ab
    print 'R_ab_rev =', R_ab_rev

    e__ab_x = R_ab * ex * R_ab_rev
    e__ab_y = R_ab * ey * R_ab_rev
    e__ab_z = R_ab * ez * R_ab_rev

    print 'e_ab_x =', e__ab_x
    print 'e_ab_y =', e__ab_y
    print 'e_ab_z =', e__ab_z

    R_phi = cos(phi/2)-(ex^ey)*sin(phi/2)
    R_phi_rev = R_phi.rev()

    print R_phi
    print R_phi_rev

    e_phi = (R_phi * ey * R_phi.rev())

    print e_phi

    R_th = cos(th/2)+I*e_phi*sin(th/2)
    R_th_rev = R_th.rev()
    print R_th
    print R_th_rev

    e_r = (R_th*R_phi*ex*R_phi_rev*R_th_rev).trigsimp()

    e_th = (R_th*R_phi*ez*R_phi_rev*R_th_rev).trigsimp()

    e_phi = e_phi.trigsimp()

    print 'e_r =', e_r
    print 'e_th =', e_th
    print 'e_phi =', e_phi

    return

if __name__ == "__main__":
    main()
