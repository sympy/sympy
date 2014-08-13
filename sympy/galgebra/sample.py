import sys, os
from sympy import Symbol, cos, sin
from sympy.galgebra.ga import Ga
from sympy.galgebra.mv import Mv
import math

def init_conformal_basis():
    basis = 'no e_1 e_2 e_3 ni'
    metric = '0 0 0 0 -1, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, -1 0 0 0 0'
    (cf3d, no, e1, e2, e3, ni) = Ga.build(basis, g=metric)
    return (no, e1, e2, e3, ni)


def make_conformal_point(epnt, basis):
    (no, e1, e2, e3, ni) = basis
    #return no + epnt[0]*e1 + epnt[1]*e2 + epnt[2]*e3 \
    #        + 0.5*math.sqrt(epnt[0]**2 + epnt[1]**2 + epnt[2]**2)*ni
    return no + epnt[0]*e1 + epnt[1]*e2 + epnt[2]*e3 \
            + 0.5*(epnt[0]**2 + epnt[1]**2 + epnt[2]**2)*ni


def gmp(q,a):
    # qAq-1 = a'
    a_p = q*a*q.inv()
    print ''
    print "Calculating a_p' = q * a * inv(q)"
    print 'q:   ', q
    print 'a:   ', a
    print 'a_p: ', a_p
    print ''
    return a_p


if __name__ == "__main__":
    """ For testing use. Typically functions will be called separately. """
    # Ideally something could be set here to lower the precision in the GA calculations.
    b = init_conformal_basis()
    (no, e1, e2, e3, ni) = b
    p1a = make_conformal_point([-8649.51, 4688.51, 600.0], b)
    p1b = make_conformal_point([4557.58, 679.734, 302.5], b)
    p2a = make_conformal_point([-8625.71, 4720.65, 600.0], b)
    p2b = make_conformal_point([4545.3, 641.665, 302.5], b)

    v1 = 1.0 + 0.0*(no+e1+e2+e3+ni) # default

    # First move - works
    marker_1_gmp = gmp(v1,p1a)
    marker_2_gmp = gmp(v1,p1b)
    #(marker_1_gmp - marker_2_gmp)^ni != marker_1_gmp^ni - marker_2_gmp^ni  !!!!operator precedence!!!!
    vec_diff = (marker_1_gmp - marker_2_gmp)^ni
    new_m2_gmp = marker_2_gmp + vec_diff
    print ''
    print 'Successful translation'
    
    # translator = (1 - (d/2)*ni)
    #v2 = (1 - 0.5*vec_diff*ni) * v1  <-- old line with typo
    v2 =  (1 - 0.5*vec_diff) * v1 # translate component versor
    print ''
    print 'v2:  ', v2
    print ''

    # Second move
    to_pt = gmp(v1,p2a)
    from_pt = gmp(v2,p2b)  # THIS PRODUCES THE LINE I REFERRED TO IN POST where
                           # latter components could be treated as zero if less demanding tolerance was set.
    center = marker_1_gmp
    to_diff = to_pt - center
    from_diff = from_pt - center
    print ''
    print "To_diff: ", to_diff
    print "From_diff: ", from_diff
    print ''
    print 'To_Diff_Norm:   ', to_diff.norm()
    print 'From_Diff_Norm: ', from_diff.norm()  # This fails because of the 1E-9










