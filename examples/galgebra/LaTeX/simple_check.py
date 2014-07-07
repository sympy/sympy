from sympy.galgebra.printer import xpdf, Get_Program, Print_Function, Format
from sympy.galgebra.ga import Ga

def basic_multivector_operations_3D():
    Print_Function()

    (g3d,ex,ey,ez) = Ga.build('e*x|y|z')

    print 'g_{ij} =',g3d.g

    A = g3d.mv('A','mv')

    A.Fmt(1,'A')
    A.Fmt(2,'A')
    A.Fmt(3,'A')

    A.even().Fmt(1,'%A_{+}')
    A.odd().Fmt(1,'%A_{-}')

    X = g3d.mv('X','vector')
    Y = g3d.mv('Y','vector')

    X.Fmt(1,'X')
    Y.Fmt(1,'Y')

    (X*Y).Fmt(2,'X*Y')
    (X^Y).Fmt(2,'X^Y')
    (X|Y).Fmt(2,'X|Y')
    return

def basic_multivector_operations_2D():
    Print_Function()

    (g2d,ex,ey) = Ga.build('e*x|y')

    print 'g_{ij} =',g2d.g

    X = g2d.mv('X','vector')
    A = g2d.mv('A','spinor')

    X.Fmt(1,'X')
    A.Fmt(1,'A')

    (X|A).Fmt(2,'X|A')
    (X<A).Fmt(2,'X<A')
    (A>X).Fmt(2,'A>X')
    return

def dummy():
    return

def main():
    Get_Program()
    Format()

    basic_multivector_operations_3D()
    basic_multivector_operations_2D()

    xpdf()
    return

if __name__ == "__main__":
    main()
