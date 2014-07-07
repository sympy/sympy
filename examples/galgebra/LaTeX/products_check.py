from sympy import symbols
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Format,xpdf

def main():
    #Format()

    coords = (x,y,z) = symbols('x y z',real=True)

    (o3d,ex,ey,ez) = Ga.build('e*x|y|z',g=[1,1,1],coords=coords)

    s = o3d.mv('s','scalar')
    v = o3d.mv('v','vector')
    b = o3d.mv('b','bivector')

    print r'#3D Orthogonal Metric\newline'

    print '#Multvectors:'
    print 's =',s
    print 'v =',v
    print 'b =',b

    print '#Products:'

    X = ((s,'s'),(v,'v'),(b,'b'))

    for xi in X:
        print ''
        for yi in X:
            print xi[1]+'*'+yi[1]+' =',xi[0]*yi[0]
            print xi[1]+'^'+yi[1]+' =',xi[0]^yi[0]
            if xi[1] != 's' and yi[1] != 's':
                print xi[1]+'|'+yi[1]+' =',xi[0]|yi[0]
            print xi[1]+'<'+yi[1]+' =',xi[0]<yi[0]
            print xi[1]+'>'+yi[1]+' =',xi[0]>yi[0]

    fs = o3d.mv('s','scalar',f=True)
    fv = o3d.mv('v','vector',f=True)
    fb = o3d.mv('b','bivector',f=True)

    print '#Multivector Functions:'

    print 's(X) =',fs
    print 'v(X) =',fv
    print 'b(X) =',fb

    print '#Products:'

    fX = ((o3d.grad,'grad'),(fs,'s'),(fv,'v'),(fb,'b'))

    for xi in fX:
        print ''
        for yi in fX:
            if xi[1] == 'grad' and yi[1] == 'grad':
                pass
            else:
                print xi[1]+'*'+yi[1]+' =',xi[0]*yi[0]
                print xi[1]+'^'+yi[1]+' =',xi[0]^yi[0]
                if xi[1] != 's' and yi[1] != 's':
                    print xi[1]+'|'+yi[1]+' =',xi[0]|yi[0]
                print xi[1]+'<'+yi[1]+' =',xi[0]<yi[0]
                print xi[1]+'>'+yi[1]+' =',xi[0]>yi[0]


    (g2d,ex,ey) = Ga.build('e',coords=(x,y))

    print r'#General 2D Metric\newline'
    print '#Multivector Functions:'

    s = g2d.mv('s','scalar',f=True)
    v = g2d.mv('v','vector',f=True)
    b = g2d.mv('v','bivector',f=True)

    print 's(X) =',s
    print 'v(X) =',v
    print 'b(X) =',b

    X = ((g2d.grad,'grad'),(s,'s'),(v,'v'))

    H =  v*g2d.grad
    print H.terms

    """
    print '#Products:'

    for xi in X:
        print ''
        for yi in X:
            if xi[1] == 'grad' and yi[1] == 'grad':
                pass
            else:
                print xi[1]+'*'+yi[1]+' =',xi[0]*yi[0]
                print xi[1]+'^'+yi[1]+' =',xi[0]^yi[0]
                if xi[1] != 's' and yi[1] != 's':
                    print xi[1]+'|'+yi[1]+' =',xi[0]|yi[0]
                print xi[1]+'<'+yi[1]+' =',xi[0]<yi[0]
                print xi[1]+'>'+yi[1]+' =',xi[0]>yi[0]
    """
    #xpdf(paper='landscape')
    return

if __name__ == "__main__":
    main()
