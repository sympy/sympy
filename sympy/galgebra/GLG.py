#Lie Algebras
from sympy import symbols
from sympy.galgebra.ga import Ga
from sympy.galgebra.mv import Com
from sympy.galgebra.printer import Format, xpdf

Estr1 = r'E_{ij}'
Fstr1 = r'F_{ij}'
Kstr1= r'K_{i}'
Estr2 = r'E_{mn}'
Fstr2 = r'F_{mn}'
Kstr2= r'K_{m}'
lm = '% '
sp = r'\:\:'
eq = ' ='
x = r' \times '

def E(i,j):
    global e,eb
    B = e[i]*e[j] - eb[i]*eb[j]
    Bstr = 'E_{' + i + j + '} '
    print '%' + Bstr + ' =',B
    return Bstr, B

def F(i,j):
    global e,eb
    B = e[i]*eb[j] - eb[i]*e[j]
    Bstr = 'F_{' + i + j + '} '
    print '%' + Bstr + ' =',B
    return Bstr, B

def K(i):
    global e,eb
    B = e[i]*eb[i]
    Bstr = 'K_{' + i + '} '
    print '%' + Bstr + ' =',B
    return Bstr, B

def ComP(A,B):
    AxB = Com(A[1],B[1])
    AxBstr = '%' + A[0] + r' \times ' + B[0] + ' ='
    print AxBstr, AxB
    return

Format()

(glg,ei,ej,em,en,ebi,ebj,ebm,ebn) = Ga.build(r'e_i e_j e_k e_l \bar{e}_i \bar{e}_j \bar{e}_k \bar{e}_l',g=[1,1,1,1,-1,-1,-1,-1])

e = {'i':ei,'j':ej,'k':em,'l':en}
eb = {'i':ebi,'j':ebj,'k':ebm,'l':ebn}

print r'#\center{General Linear Group of Order $n$\newline}'
print r'#Lie Algebra Generators: $1\le i < j \le n$ and $1 \le i < l \le n$'

Eij = E('i','j')
Fij = F('i','j')
Ki = K('i')
Eil = E('i','l')
Fil = F('i','l')

print r'#Non Zero Commutators'

ComP(Eij,Fij)
ComP(Eij,Ki)
ComP(Fij,Ki)
ComP(Eij,Eil)
ComP(Fij,Fil)
ComP(Fij,Eil)

xpdf(paper='letter',pt='12pt',debug=True,prog=True)


