from sympy import symbols, sin, cos
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Format, xpdf, Eprint, Print_Function, Get_Program

def main():
    Print_Function()

    (a, b, c) = abc = symbols('a,b,c',real=True)
    (o3d, ea, eb, ec) = Ga.build('e_a e_b e_c', g=[1, 1, 1], coords=abc)
    grad = o3d.grad

    x = symbols('x',real=True)
    A = o3d.lt([[x*a*c**2,x**2*a*b*c,x**2*a**3*b**5],\
                [x**3*a**2*b*c,x**4*a*b**2*c**5,5*x**4*a*b**2*c],\
                [x**4*a*b**2*c**4,4*x**4*a*b**2*c**2,4*x**4*a**5*b**2*c]])
    print 'A =',A

    v = a*ea+b*eb+c*ec

    print 'v =',v

    f = v|A(v)

    print r'%f = v\cdot \f{A}{v} =',f

    (grad * f).Fmt(3,r'%\nabla f')

    Av = A(v)

    print r'%\f{A}{v} =', Av

    (grad * Av).Fmt(3,r'%\nabla \f{A}{v}')

    return

def dummy():
    return

if __name__ == "__main__":
    #Eprint()
    Format()
    Get_Program()
    main()
    xpdf()
