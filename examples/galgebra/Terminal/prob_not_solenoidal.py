"""
import sys
sys.path.append('/home/brombo/MathPhysics/sympy/sympy/galgebra')
"""
from sympy import symbols,sqrt
from sympy.galgebra.printer import Eprint
from sympy.galgebra.ga import Ga

def main():
    Eprint()

    X = (x,y,z) = symbols('x y z',real=True)
    (o3d,ex,ey,ez) = Ga.build('e_x e_y e_z',g=[1,1,1],coords=(x,y,z))

    A = x*(ey^ez) + y*(ez^ex) + z*(ex^ey)
    print 'A =', A
    print 'grad^A =',(o3d.grad^A).simplify()
    print

    f = o3d.mv(1/sqrt(x**2 + y**2 + z**2))
    print 'f =', f
    print 'grad*f =',(o3d.grad*f).simplify()
    print

    B = f*A
    print 'B =', B
    print

    Curl_B = o3d.grad^B

    print 'grad^B =', Curl_B.simplify()

    return

if __name__ == "__main__":
    main()
