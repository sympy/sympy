#Lie Algebras
from sympy import symbols
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import GaLatexPrinter
from sympy.galgebra.printer import Format, xpdf

class LieAlgebra(object):

    def __init__(self,n):
        self.n = n
        e = ''
        ebar = ''
        for i in range(1,n+1):
            e += ' e_' + str(i)
            if GaLatexPrinter.latex_flg:
                ebar += r' \bar{e}_' + str(i)
            else:
                ebar += ' ebar_' + str(i)

        g = n * [1] + n * [-1]
        basis = e[1:] + ebar

        self.Ga = Ga(basis, g=g)
        self.basis = self.Ga.mv()

        self.e = self.basis[:n]
        self.ebar = self.basis[n:]

        self.w = []
        self.wstar = []

        for i in range(n):
            self.w.append(self.e[i] + self.ebar[i])
            self.wstar.append(self.e[i] - self.ebar[i])

        self.Nu_bais = self.w + self.wstar
        self.Eij = []
        self.Fij = []
        self.Ki = []
        for i in range(n):
            self.Ki.append(self.e[i] * self.ebar[i])
            print r'%F_{'+str(i)+'} =',self.Ki[-1] * self.Ki[-1].rev()
            for j in range(i):
                self.Eij.append(self.e[i] * self.e[j] - self.ebar[i] * self.ebar[j])
                self.Fij.append(self.e[i] * self.ebar[j] - self.ebar[i] * self.e[j])
                print r'%E_{'+str(i)+str(j)+'} =',self.Eij[-1] * self.Eij[-1].rev()
                print r'%F_{'+str(i)+str(j)+'} =',self.Fij[-1] * self.Fij[-1].rev()

        print 'K_{i} =',self.Ki
        print 'E_{ij} =',self.Eij
        print 'F_{ij} =',self.Fij

        E = self.Eij[0]/2

        for i in range(2*n):
            print E
            E *= self.Eij[0]/2



if __name__ == "__main__":

    Format()
    la = LieAlgebra(6)
    print la.e
    print la.ebar
    xpdf()
