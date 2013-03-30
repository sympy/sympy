import sys

from sympy import symbols,sin,cos
from sympy.ga.ga import MV,Format
from sympy.ga.ga_print import xdvi,Get_Program,Print_Function

Format()

def derivatives_in_spherical_coordinates():
    Print_Function()
    X = (r,th,phi) = symbols('r theta phi')
    curv = [[r*cos(phi)*sin(th),r*sin(phi)*sin(th),r*cos(th)],[1,r,r*sin(th)]]
    (er,eth,ephi,grad) = MV.setup('e_r e_theta e_phi',metric='[1,1,1]',coords=X,curv=curv)

    f = MV('f','scalar',fct=True)
    A = MV('A','vector',fct=True)
    B = MV('B','grade2',fct=True)

    print 'f =',f
    print 'A =',A
    print 'B =',B

    print 'grad*f =',grad*f
    print 'grad|A =',grad|A
    print '-I*(grad^A) =',-MV.I*(grad^A)
    print 'grad^B =',grad^B
    return
def dummy():
    return


Get_Program(True)
derivatives_in_spherical_coordinates()
xdvi()
