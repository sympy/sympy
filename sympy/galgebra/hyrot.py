from sympy import symbols, sinh, cosh
from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga

def Com(x,y):
    return x * y - y * x

def LieBasis(n):
    g = n * [1] + n * [-1]
    basis = ''
    for i in range(1,n+1):
        basis += 'e_' +str(i) + ' '
    for i in range(1,n+1):
        basis += r'\bar{e}_' +str(i) + ' '
    basis = basis[:-1]

    LieGA = Ga(basis,g=g)
    bases = LieGA.mv()
    e = bases[:n]
    ebar = bases[n:]
    print e
    print ebar

    print LieGA.g

    E = []
    F = []
    K = []
    indexes = []

    for i in range(n):
        K.append(e[i]*ebar[i])
        for j in range(n):
            if i < j:
                indexes.append((i+1,j+1))
                E.append(e[i]*e[j]-ebar[i]*ebar[j])
                F.append(e[i]*ebar[j]-ebar[i]*e[j])

    print indexes
    print 'E =',E
    print 'F =',F
    print 'K =',K,'\n'

    for k in range(len(indexes)):
        k_i = indexes[k][0]
        k_j = indexes[k][1]
        for l in range(len(indexes)):
            l_i = indexes[l][0]
            l_j = indexes[l][1]
            print 'E_'+str(k_i)+str(k_j)+' x F_'+str(l_i)+str(l_j)+' = '+str(Com(E[k],F[l]))

    return


LieBasis(6)

